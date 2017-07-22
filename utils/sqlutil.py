from __future__ import print_function

# This file is part of astrolibpy
#
#	 astrolibpy is free software: you can redistribute it and/or modify
#	 it under the terms of the GNU General Public License as published by
#	 the Free Software Foundation, either version 3 of the License, or
#	 (at your option) any later version.
#
#	astrolibpy is distributed in the hope that it will be useful,
#	 but WITHOUT ANY WARRANTY; without even the implied warranty of
#	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	 GNU General Public License for more details.
#
#	 You should have received a copy of the GNU General Public License
#	 along with astrolibpy.	 If not, see <http://www.gnu.org/licenses/>.


import types
import numpy
import sys
import numpy as np
import time
import psycopg2
import threading
try:
	from cStringIO import StringIO
except ImportError:
	# for python3
	from io import BytesIO as  StringIO 
try:
	import queue
except:
	import Queue as queue
import collections
import warnings
try:
	dictclass = collections.OrderedDict
except AttributeError:
	try:
		import ordereddict
		dictclass = ordereddict.OrderedDict
	except ImportError:
		warnings.warn('No ordered dict library found')
		dictclass = dict
from numpy.core import numeric as sb
from numpy.core import numerictypes as nt
from select import select
from psycopg2.extensions import POLL_OK, POLL_READ, POLL_WRITE

def wait_select_inter(conn):
	# Make the queries interruptable by Ctrl-C
	# Taken from http://initd.org/psycopg/articles/2014/07/20/cancelling-postgresql-statements-python/
	while True:
		try:
			state = conn.poll()
			if state == POLL_OK:
				break
			elif state == POLL_READ:
				select([conn.fileno()], [], [])
			elif state == POLL_WRITE:
				select([], [conn.fileno()], [])
			else:
				raise conn.OperationalError(
					"bad state from poll: %s" % state)
		except KeyboardInterrupt:
			conn.cancel()
			# the loop will be broken by a server error
			continue
psycopg2.extensions.set_wait_callback(wait_select_inter)

def getConnection(db=None, driver=None, user=None,
				  password=None, host=None, port=5432, timeout=None):
	if driver == 'psycopg2':
		conn_str = "dbname=%s host=%s port=%d" % (db, host, port)
		if user is not None:
			conn_str = conn_str + ' user=%s' % user
		if password is not None:
			conn_str = conn_str + ' password=%s' % password
		conn = psycopg2.connect(conn_str)
	elif driver == 'sqlite3':
		import sqlite3
		if timeout is None:
			timeout = 5
		conn = sqlite3.connect(db, timeout=timeout)
	else:
		raise Exception("Unknown driver")
	return conn


def getCursor(conn, driver=None, preamb=None, notNamed=False):
	if driver == 'psycopg2':
		cur = conn.cursor()
		if preamb is not None:
			cur.execute(preamb)
		else:
			cur.execute('set cursor_tuple_fraction TO 1')
			# this is required because otherwise PG may decide to execute a
			# different plan
		if notNamed:
			return cur
		cur = conn.cursor(name='sqlutilcursor')
		cur.arraysize = 100000
	elif driver == 'sqlite3':
		cur = conn.cursor()
	return cur


def fromrecords(recList, dtype=None, intNullVal=None):
	""" This function was taken from np.core.records and updated to
			support conversion null integers to intNullVal
	"""

	nfields = len(recList[0])
	shape = None 
	descr = sb.dtype((np.core.records.record, dtype))

	try:
		retval = sb.array(recList, dtype=descr)
	except TypeError:  # list of lists instead of list of tuples
		shape = (len(recList),)
		_array = np.core.records.recarray(shape, descr)
		try:
			for k in range(_array.size):
				_array[k] = tuple(recList[k])
		except TypeError:
			convs = []
			ncols = len(dtype.fields)
			for _k in dtype.names:
				_v = dtype.fields[_k]
				if _v[0] in [np.int16, np.int32, np.int64]:
					convs.append(lambda x: intNullVal if x is None else x)
				else:
					convs.append(lambda x: x)
			convs = tuple(convs)
			convF = lambda x: [convs[_](x[_]) for _ in range(ncols)]

			for k in range(k, _array.size):
				try:
					_array[k] = tuple(recList[k])
				except TypeError:
					_array[k] = tuple(convF(recList[k]))
		return _array
	else:
		if shape is not None and retval.shape != shape:
			retval.shape = shape

	res = retval.view(numpy.core.records.recarray)

	return res


def __converter(qIn, qOut, endEvent, dtype, intNullVal):
	""" Convert the input stream of tuples into numpy arrays """
	while(not endEvent.is_set()):
		try:
			tups = qIn.get(True, 0.1)
		except queue.Empty:
			continue
		try:
			res = fromrecords(tups, dtype=dtype, intNullVal=intNullVal)
		except:
			print ('Failed to convert input data into array')
			endEvent.set()
			raise
		qOut.put(res)


def get(query, params=None, db="wsdb", driver="psycopg2", user=None,
		password=None, host='localhost', preamb=None,
		conn=None, port=5432,
		strLength=10, timeout=None, notNamed=False,
		asDict=False, intNullVal=-9999):
	'''Executes the sql query and returns the tuple or dictionary with the numpy arrays.

	Parameters:
	------
	query : string with the query you want to execute, can include question 
			marks to refer to query parameters
	params : tuple with query parameters
	conn : the connection object to the DB (optional) to avoid reconnecting
	asDict : boolean to retrieve the results as a dictionary with column names
			as keys
	strLength : all the strings will be truncated to this length
	intNullVal : all the integer columns with nulls will have null replaced by
				 this value
	db : string with the name of the database
	driver : the sql driver to be used (psycopg2 and sqlite3 are supported)
	user : user name for the DB connection
	password : DB connection password
	host : hostname of the database
	port : port of the database 
	preamb: bit of SQL code to be executed before the query

	Example:
	>>> a, b, c = sqlutil.get('select ra,dec,d25 from rc3')
	You can also use the parameters in your query:
	Example:
	>>> a, b = squlil.get('select ra,dec from rc3 where name=?',"NGC 3166")
	'''
	__pgTypeHash = {
		16: bool,
		18: str,
		20: 'i8',
		21: 'i2',
		23: 'i4',
		1007: 'i4',
		25: '|S%d' % strLength,
		700: 'f4',
		701: 'f8',
			1042: '|S%d' % strLength,  # character()
			1043: '|S%d' % strLength,  # varchar
			1700: 'f8',	 # numeric
			1114: '<M8[us]', # timestamp
			1082: '<M8[us]'	 # date
	}

	connSupplied = (conn is not None)
	if not connSupplied:
		conn = getConnection(db=db, driver=driver, user=user, password=password,
							 host=host, port=port, timeout=timeout)
	try:
		cur = getCursor(conn, driver=driver, preamb=preamb, notNamed=notNamed)

		if params is None:
			res = cur.execute(query)
		else:
			res = cur.execute(query, params)

		qIn = queue.Queue(1)
		qOut = queue.Queue()
		endEvent = threading.Event()
		nrec = 0  # keeps the number of arrays sent to the other thread
		# minus number received
		reslist = []
		proc = None
		colNames = []
		if driver == 'psycopg2':
			try:
				while(True):
					tups = cur.fetchmany()
					if nrec == 0:
						desc = cur.description
						typeCodes = [_tmp.type_code for _tmp in desc]
						colNames = [_tmp.name for _tmp in cur.description]
						dtype = numpy.dtype(
							[('a%d' % _i, __pgTypeHash[_t]) for _i, _t in enumerate(typeCodes)])
						proc = threading.Thread(target=__converter, args=(
							qIn, qOut, endEvent, dtype, intNullVal))
						proc.start()

					if tups == []:
						break
					qIn.put(tups)
					nrec += 1
					try:
						reslist.append(qOut.get(False))
						nrec -= 1
					except queue.Empty:
						pass
				try:
					while(nrec != 0):
						try:
							reslist.append(qOut.get(True, 0.1))
							nrec -= 1
						except:
							if endEvent.is_set():
								raise Exception('Child thread failed')
				except queue.Empty:
					pass
				endEvent.set()
			except BaseException:
				endEvent.set()
				if proc is not None:
					# notice that here the timeout is larger than the timeout
					proc.join(0.2)
					# in the converter process
					if proc.is_alive():
						proc.terminate()
				raise
			proc.join()
			if reslist == []:
				nCols = len(desc)
				res = numpy.array([],
								  dtype=numpy.dtype(
									  [('a%d' % i, 'f') for i in range(nCols)])
								  )
			else:
				res = numpy.concatenate(reslist)

		elif driver == 'sqlite3':
			tups = cur.fetchall()
			if len(tups) > 0:
				_cast = {types.BooleanType: numpy.bool,
						 types.IntType: numpy.int32,
						 types.LongType: numpy.int64,
						 types.FloatType: numpy.float64,
						 types.StringType: numpy.str,
						 types.UnicodeType: numpy.str}
				try:
					typelist = [_cast[type(tmp)] for tmp in tups[0]]
				except KeyError:
					raise Exception("Unknown datatype")
				res = numpy.core.records.array(tups)
			else:
				return None

		res = [res[tmp] for tmp in res.dtype.names]

	except BaseException:
		try:
			conn.rollback()
		except:
			pass
		if not connSupplied:
			try:
				conn.close()  # do not close if we were given the connection
			except:
				pass
		raise

	cur.close()

	if asDict:
		resDict = dictclass()
		repeats = {}
		for _n, _v in zip(colNames, res):
			if _n in resDict:
				curn = _n + '_%d'%(repeats[_n])
				repeats[_n]+=1
				warnings.warn(('Column name %s is repeated in the output, '+
								'new name %s assigned')%(_n,curn))
			else:
				repeats[_n]=1
				curn = _n
			resDict[curn] = _v
		res = resDict
	return res

def execute(query, params=None, db="wsdb", driver="psycopg2", user=None,
			password=None, host='locahost',
			conn=None, preamb=None, timeout=None,
			noCommit=False):
	"""Execute a given SQL command without returning the results"""
	connSupplied = (conn is not None)
	if not connSupplied:
		conn = getConnection(db=db, driver=driver, user=user, password=password,
							 host=host, timeout=timeout)
	try:
		cur = getCursor(conn, driver=driver, preamb=preamb, notNamed=True)
		if params is not None:
			cur.execute(query, params)
		else:
			# sqlite3 doesn't like params here...
			cur.execute(query)
	except BaseException:
		try:
			conn.rollback()
		except Exception:
			pass
		if not connSupplied:
			try:
				conn.close()  # do not close if we were given the connection
			except:
				pass
		raise
	cur.close()
	if not noCommit:
		conn.commit()
	if not connSupplied:
		conn.close()  # do not close if we were given the connection


def __create_schema(tableName, arrays, names, temp=False):
	hash = dict([
		(np.int32, 'integer'),
		(np.int64, 'bigint'),
		(np.float32, 'real'),
		(np.float64, 'double precision'),
		(np.string_, 'varchar'),
	])
	if temp:
		temp = 'temporary'
	else:
		temp = ''
	outp = 'create %s table %s ' % (temp, tableName)
	outp1 = []
	for arr, name in zip(arrays, names):
		outp1.append(name + ' ' + hash[arr.dtype.type])
	return outp + '(' + ','.join(outp1) + ')'


def __print_arrays(arrays, f):
	hash = dict([
		(np.int32, '%d'),
		(np.int64, '%d'),
		(np.float32, '%.18e'),
		(np.float64, '%.18e'),
		(np.string_, '%s')
	])
	fmt = [hash[x.dtype.type] for x in arrays]
	recarr = np.rec.fromarrays(arrays)
	np.savetxt(f, recarr, fmt=fmt)


def upload(tableName, arrays, names, db="wsdb", driver="psycopg2", user=None,
		   password=None, host='locahost',
		   conn=None, preamb=None, timeout=None,
		   noCommit=False, temp=False,
		   analyze=False, createTable=True):
	""" Upload the data stored in the tuple of arrays in the DB

	Example:
	>>> x = np.arange(10)
	>>> y = x**.5
	>>> sqlutil.upload('mytable',(x,y),('xcol','ycol'))
	"""
	arrays = [np.asarray(_) for _ in arrays]
	connSupplied = (conn is not None)
	if not connSupplied:
		conn = getConnection(db=db, driver=driver, user=user, password=password,
							 host=host, timeout=timeout)
	try:
		cur = getCursor(conn, driver=driver, preamb=preamb, notNamed=True)
		if createTable:
			query1 = __create_schema(tableName, arrays, names, temp=temp)
			cur.execute(query1)
		f = StringIO()
		__print_arrays(arrays, f)
		f.seek(0)
		try:
			thread = psycopg2.extensions.get_wait_callback()
			psycopg2.extensions.set_wait_callback(None)
			cur.copy_from(f, tableName, sep=' ', columns=names)
		finally:
			psycopg2.extensions.set_wait_callback(thread)
	except BaseException:
		try:
			conn.rollback()
		except Exception:
			pass
		if not connSupplied:
			try:
				conn.close()  # do not close if we were given the connection
			except:
				pass
		raise
	if analyze:
		cur.execute('analyze %s' % tableName)
	cur.close()
	if not noCommit:
		conn.commit()
	if not connSupplied:
		conn.close()  # do not close if we were given the connection


def local_join(query, tableName, arrays, names, db="wsdb", driver="psycopg2", user=None,
			   password=None, host='locahost',
			   port=5432,
			   conn=None, preamb=None, timeout=None,
			   strLength=20, asDict=False):
	""" Join the data from python with the data in the database
	This command first uploads the data in the DB and then runs a 
	user specified query.

	Parameters
	----------
	query : String with the query to be executed 
	tableName : The name of the temporary table that is going to be created
	arrays : The tuple with list of arrays with the data to be loaded in the DB
	names : The tuple with the column names for the user table

	Example: 
	>>> x = np.arange(10)
	>>> y = x**.5
	>>> sqlutil.local_join('select * from mytable as m, sometable as s where s.id=m.xcol', 
													'mytable',(x,y),('xcol','ycol'))
	"""

	connSupplied = (conn is not None)
	if not connSupplied:
		conn = getConnection(db=db, driver=driver, user=user, password=password,
							 host=host, timeout=timeout, port=port)

	upload(tableName, arrays, names, conn=conn, noCommit=True, temp=True,
		   analyze=True)
	res = get(query, conn=conn, preamb=preamb, strLength=strLength,
			  asDict=asDict)
	conn.rollback()

	if not connSupplied:
		conn.close()
	return res


# This file is part of astrolibpy
#
#    astrolibpy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#   astrolibpy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with astrolibpy.  If not, see <http://www.gnu.org/licenses/>.


import types
import numpy, sys, numpy as np
import time,psycopg2
import threading, Queue
import collections, warnings
try:
	dictclass = collections.OrderedDict
except AttributeError:
	try:
		import ordereddict
		dictclass = ordereddict.OrderedDict
	except ImportError:
		warnings.warn('No ordered dict library found')
		dictclass = dict

def getConnection( db=None, driver=None, user=None,
						password=None, host=None,port=5432, timeout=None):
	if driver=='psycopg2':
		import psycopg2
		conn_str = "dbname=%s host=%s port=%d"%(db,host,port)
		if user is not None:
			conn_str = conn_str+ ' user=%s'%user
		if password is not None:
			conn_str = conn_str+ ' password=%s'%password
		conn = psycopg2.connect(conn_str)
	elif driver=='sqlite3':
		import sqlite3
		if timeout is None:
			timeout = 5
		conn = sqlite3.connect(db, timeout=timeout)
	else: 
		raise Exception("Unknown driver")
	return conn					

def getCursor(conn, driver=None, preamb=None, notNamed=False):
	if driver=='psycopg2':
		cur = conn.cursor()
		if preamb is not None:
			cur.execute(preamb)
		else:
			cur.execute('set cursor_tuple_fraction TO 1') 
			# this is required because otherwise PG may decide to execute a different plan
		if notNamed:
			return cur
		cur = conn.cursor(name='sqlutilcursor')
		cur.arraysize=100000
	elif driver=='sqlite3':
		cur = conn.cursor()
	return cur

def __converter(qIn, qOut, endEvent, dtype):
	while(not endEvent.is_set()):
		try:
			tups = qIn.get(True,0.1)
		except Queue.Empty:
			continue
		try:
			res=numpy.core.records.array(tups,dtype=dtype)
		except:
			print 'Failed to convert input data into array'
			endEvent.set()
			raise
		qOut.put(res)


def get(query, params=None, db="wsdb", driver="psycopg2", user=None,
						password=None, host='localhost', preamb=None,
						getConn=False, conn=None, maskNull=False, port=5432,
						strLength=10, timeout=None, notNamed=False, asDict=False):
	'''This program executes the sql query and returns 
	the tuple of the numpy arrays.
	Example:
	a,b, c = sqlutil.get('select ra,dec,d25 from rc3')
	You can also use the parameters in your query:
	Example:
	a,b = squlil.get('select ra,dec from rc3 where name=?',"NGC 3166")
	'''
	__pgTypeHash = {16:bool,18:str,20:'i8',21:'i2',23:'i4',25:'|S%d'%strLength,700:'f4',701:'f8',
		1043:'|S%d'%strLength,#varchar
		1700:'f8' #numeric
		} 

	connSupplied = (conn is not None)
	if not connSupplied:
		conn = getConnection(db=db, driver=driver, user=user, password=password,
				host=host, port=port, timeout=timeout)
	try:
		cur = getCursor(conn, driver=driver, preamb=preamb,notNamed=notNamed)

		if params is None:
			res = cur.execute(query)
		else:
			res = cur.execute(query, params)
		
		qIn = Queue.Queue(1)
		qOut = Queue.Queue()
		endEvent = threading.Event()
		nrec = 0 ## keeps the number of arrays sent to the other thread
				##  minus number received
		reslist=[]
		proc = None
		colNames = []
		if driver=='psycopg2':
			try:
				while(True):
					tups = cur.fetchmany()
					if nrec==0:
						desc = cur.description
						typeCodes = [_tmp.type_code for _tmp in desc]
						colNames = [_tmp.name for _tmp in cur.description]
						dtype=numpy.dtype([('a%d'%_i,__pgTypeHash[_t]) for _i,_t in enumerate(typeCodes)])				
						proc = threading.Thread(target=__converter, args = (qIn, qOut, endEvent,dtype))
						proc.start()

					if tups == []:
						break
					qIn.put(tups)
					nrec += 1
					try:
						reslist.append(qOut.get(False))
						nrec -= 1
					except Queue.Empty:
						pass
				try:
					while(nrec!=0):
						try:
							reslist.append(qOut.get(True, 0.1))
							nrec-=1
						except:
							if endEvent.is_set():
								raise Exception('Child thread failed')
				except Queue.Empty:
					pass
				endEvent.set()
			except BaseException:
				endEvent.set()
				if proc is not None:
					proc.join(0.2) # notice that here the timeout is larger than the timeout
									# in the converter process
					if proc.is_alive():
						proc.terminate()
				raise
			proc.join()
			if reslist == []:
				nCols = len(desc)
				res = numpy.array([],
					dtype=numpy.dtype([('a%d'%i,'f') for i in range(nCols)])
									)				
			else:
				res = numpy.concatenate(reslist)

		elif driver=='sqlite3':
			tups=cur.fetchall()
			if len(tups)>0:
				_cast = {types.BooleanType: numpy.bool,
					types.IntType: numpy.int32,
					types.LongType: numpy.int64,
					types.FloatType: numpy.float64,
					types.StringType: numpy.str,
					types.UnicodeType: numpy.str}
				try:
					typelist=[_cast[type(tmp)] for tmp in tups[0]]
				except KeyError:
					raise Exception("Unknown datatype")
				res = numpy.core.records.array(tups)
			else:
				return None

		res=[res[tmp] for tmp in res.dtype.names]
		if maskNull:
			for i in range(len(res)):
				if res[i].dtype==numpy.object:
					res[i]=res[i].astype(numpy.float)

	except BaseException:
		try:
			conn.rollback()
		except Exception:
			pass
		if not connSupplied:
			try:
				conn.close() # do not close if we were given the connection
			except:
				pass
		raise

	cur.close()
	#conn.commit()
	if asDict:
		resDict = dictclass()
		for _n, _v in zip(colNames,res):
			resDict[_n]=_v
		res=resDict	
	if not getConn:
		if not connSupplied:

			conn.close() # do not close if we were given the connection
		return res
	else:
		return conn,res

def execute(query, params=None, db="wsdb", driver="psycopg2", user=None,
										password=None, host='locahost',
										conn=None, preamb=None, timeout=None,
										noCommit=False):
	"Execute a given SQL command without returning the results"
	connSupplied = (conn is not None)
	if not connSupplied:
		conn = getConnection(db=db,driver=driver,user=user,password=password,
				host=host, timeout=timeout)
	try:
		cur = getCursor(conn, driver=driver, preamb=preamb, notNamed=True)
		
		cur.execute(query, params)
	except BaseException:
		try:
			conn.rollback()
		except Exception:
			pass
		if not connSupplied:
			try:
				conn.close() # do not close if we were given the connection
			except:
				pass
		raise
	cur.close()
	if not noCommit:
		conn.commit()
	if not connSupplied:
		conn.close() # do not close if we were given the connection


def create_schema(tableName, arrays, names, temp=False):
	hash = dict([
		(np.int32,'integer'),
		(np.int64,'bigint'),
		(np.float32,'real'),
		(np.float64,'double precision'),
		(np.string_,'varchar'),		
			])
	if temp:
		temp='temporary'
	else:
		temp=''
	outp = 'create %s table %s '%(temp,tableName)
	outp1 = []
	for arr,name in zip(arrays,names):
		outp1.append(name+' '+hash[arr.dtype.type])
	return outp + '(' + ','.join(outp1)+')'

def print_arrays(arrays, f):
	hash = dict([
		(np.int32,'%d'),
		(np.int64,'%d'),
		(np.float32,'%.18e'),
		(np.float64,'%.18e'),
		(np.string_,'%s')
			])
	fmt = [ hash[x.dtype.type] for x in arrays]
	recarr = np.rec.fromarrays(arrays)	
	np.savetxt(f, recarr, fmt=fmt)
		
def upload(tableName, arrays, names, db="wsdb", driver="psycopg2", user=None,
										password=None, host='locahost',
										conn=None, preamb=None, timeout=None,
										noCommit=False, temp=False):
	""" Upload the data stored in the tuple of arrays in the DB
	x=np.arange(10)
	y=x**.5
	sqlutil.upload('mytable',(x,y),('xcol','ycol'))
	"""
	
	connSupplied = (conn is not None)
	if not connSupplied:
		conn = getConnection(db=db,driver=driver,user=user,password=password,
				host=host, timeout=timeout)
	try:
		cur = getCursor(conn, driver=driver, preamb=preamb, notNamed=True)
		query1 = create_schema(tableName, arrays, names, temp=temp)
		cur.execute(query1)
		import StringIO
		f = StringIO.StringIO()
		print_arrays(arrays, f)
		f.seek(0)
		cur.copy_from(f,tableName,sep=' ')
	except BaseException:
		try:
			conn.rollback()
		except Exception:
			pass
		if not connSupplied:
			try:
				conn.close() # do not close if we were given the connection
			except:
				pass
		raise
	cur.close()
	if not noCommit:
		conn.commit()
	if not connSupplied:
		conn.close() # do not close if we were given the connection
	

		
def local_join(query, tableName, arrays, names, db="wsdb", driver="psycopg2", user=None,
										password=None, host='locahost',
										conn=None, preamb=None, timeout=None):
	""" This function allows joining the data from python with the data in the DB,it
	first uploads the data in the DB and then run a user specified query:
	x=np.arange(10)
	y=x**.5
	sqlutil.local_join('select * from mytable as m, sometable as s where s.id=m.xcol', 
							'mytable',(x,y),('xcol','ycol'))
	"""

	connSupplied = (conn is not None)
	if not connSupplied:
		conn = getConnection(db=db,driver=driver,user=user,password=password,
				host=host, timeout=timeout)
	
	upload(tableName, arrays, names, conn=conn, noCommit=True, temp=True)
	res=get(query,conn=conn)
	
	if not connSupplied:
		conn.close()
	return res

					
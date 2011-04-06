# Copyright (C) 2009-2010 Sergey Koposov
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
import numpy, sys
import time,psycopg2
import threading, Queue


def getConnection( db=None, driver=None, user=None,
						password=None, host=None):
	if driver=='psycopg2':
		import psycopg2
		conn_str = "dbname=%s host=%s"%(db,host)
		if user is not None:
			conn_str = conn_str+ ' user=%s'%user
		if password is not None:
			conn_str = conn_str+ ' password=%s'%password
		conn = psycopg2.connect(conn_str)
	elif driver=='sqlite3':
		import sqlite3
		conn = sqlite3.connect(db)
		cur = conn.cursor()
	else: 
		raise Exception("Unknown driver")
	return conn					

def getCursor(conn, driver=None, preamb=None, notNamed=False):
	if driver=='psycopg2':
		cur = conn.cursor()
		if preamb is not None:
			cur.execute(preamb)
		if notNamed:
			return cur
		cur = conn.cursor(name='sqlutilcursor')
		cur.arraysize=100000
	elif driver=='sqlite3':
		cur = conn.cursor()
	return cur

def __converter(qIn, qOut, endEvent):
	while(not endEvent.is_set()):
		try:
			tups = qIn.get(True,0.1)
		except Queue.Empty:
			continue
		qOut.put(numpy.core.records.array(tups))


def get(query, params=None, db="wsdb", driver="psycopg2", user=None,
						password=None, host='localhost', preamb=None,
						getConn=False, conn=None, maskNull=False):
	'''This program executes the sql query and returns 
	the tuple of the numpy arrays.
	Example:
	a,b, c = sqlutil.get('select ra,dec,d25 from rc3')
	You can also use the parameters in your query:
	Example:
	a,b = squlil.get('select ra,dec from rc3 where name=?',"NGC 3166")
	'''
	connSupplied = (conn is not None)
	if not connSupplied:
		conn = getConnection(db=db,driver=driver,user=user,password=password,
				host=host)
	try:
		cur = getCursor(conn, driver=driver, preamb=preamb)

		if params is None:
			res = cur.execute(query)
		else:
			res = cur.execute(query, params)

		qIn = Queue.Queue(1)
		qOut = Queue.Queue()
		endEvent = threading.Event()
		nrec = 0
		reslist=[]
		if driver=='psycopg2':
			proc = threading.Thread(target=__converter, args = (qIn, qOut, endEvent))
			proc.start()
			try:
				while(True):
					tups = cur.fetchmany()
					if tups == []:
						break
					qIn.put(tups)
					nrec+=1
					try:
						reslist.append(qOut.get(False))
						nrec-=1
					except Queue.Empty:
						pass
				try:
					while(nrec!=0):
						reslist.append(qOut.get(True))
						nrec-=1
				except Queue.Empty:
					pass
				endEvent.set()
			except BaseException:
				endEvent.set()
				proc.join(0.2) # notice that here the timeout is larger than the timeout
								# in the converter process
				if proc.is_alive():
					proc.terminate()
				raise
			proc.join()
			if reslist == []:
				nCols = len(cur.description)
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
	conn.rollback()
	
	if not getConn:
		if not connSupplied:
			conn.close() # do not close if we were given the connection
		return res
	else:
		return conn,res

def execute(query, db="wsdb", driver="psycopg2", user=None,
										password=None, host='locahost',
										conn=None, preamb=None):
	connSupplied = (conn is not None)
	if not connSupplied:
		conn = getConnection(db=db,driver=driver,user=user,password=password,
				host=host)
	try:
		cur = getCursor(conn, driver=driver, preamb=preamb, notNamed=True)
		
		cur.execute(query)
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
	conn.commit()
	if not connSupplied:
		conn.close() # do not close if we were given the connection

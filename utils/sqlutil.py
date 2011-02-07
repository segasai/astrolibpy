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
import multiprocessing, Queue

def get(query, params=None, db="wsdb", driver="psycopg2", user=None,
										password=None, host=None):
	'''This program executes the sql query and returns 
	the tuple of the numpy arrays.
	Example:
	a,b, c = sqlutil.get('select ra,dec,d25 from rc3')
	You can also use the parameters in your query:
	Example:
	a,b = squlil.get('select ra,dec from rc3 where name=?',"NGC 3166")
	'''
	if driver=='psycopg2':
		import psycopg2
		con = psycopg2.connect("dbname=%s user=%s password=%s host=%s"%(db,user,password,host))
		cur = con.cursor(name='sqlutilcursor')
		cur.arraysize=100000

	elif driver=='sqlite3':
		import sqlite3
		con = sqlite3.connect(db)
		cur = con.cursor()
	else: 
		raise Exception("Unknown driver")
	if params==None:
		res = cur.execute(query)
	else:
		res = cur.execute(query, params)
	qIn = multiprocessing.Queue(1)
	qOut = multiprocessing.Queue()
	endEvent = multiprocessing.Event()
	endEvent1 = multiprocessing.Event()
	
	def converter():
		try:
			while(not endEvent.is_set()):
				try:
					tups = qIn.get(True,0.1)
				except Queue.Empty:
					continue
				qOut.put(numpy.core.records.array(tups))
		finally:
			endEvent1.set()

	reslist=[]
	proc = multiprocessing.Process(target=converter)
	if driver=='psycopg2':
		proc.start()
		try:
			while(True):
				tups = cur.fetchmany()
				if tups == []:
					endEvent.set()
					break
				qIn.put(tups)
				try:
					reslist.append(qOut.get(False))
				except Queue.Empty:
					pass
			try:
				endEvent1.wait()
				while(True):
					reslist.append(qOut.get(False))
			except Queue.Empty:
				pass
		except BaseException:
			ei = sys.exc_info()
			endEvent.set()
			proc.join(0.2) # notice that here the timeout is larger than the timeout
							# in the converter process
			if proc.is_alive():
				proc.terminate()
			raise ei[0], None, ei[2]
		proc.join()
		if reslist == []:
			return None
		res=numpy.concatenate(reslist)

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

	cur.close()
	con.commit()
	con.close()
	return res

def execute(query, db="wsdb", driver="psycopg2", user=None,
										password=None, host=None):
	if driver=='psycopg2':
		import psycopg2
		con = psycopg2.connect("dbname=%s user=%s password=%s host=%s"%(db,user,password,host))
	elif driver=='sqlite3':
		import sqlite3
		con = sqlite3.connect(db)
	else:
		raise Exception('Unrecognized driver')
	cur=con.cursor()
	cur.execute(query)
	con.commit()
	con.close()

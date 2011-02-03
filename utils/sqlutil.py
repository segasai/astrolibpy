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
import numpy
import time,psycopg2


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
	elif driver=='sqlite3':
		import sqlite3
		con = sqlite3.connect(db)
	else: 
		pass
#		raise Exception("Unknown driver")
	cur = con.cursor(name='mycursor')
	cur.arraysize=1000000
	if params==None:
		res = cur.execute(query)
	else:
		res = cur.execute(query, params)

	reslist=[]
	if driver=='psycopg2':
		while(True):
			tups = cur.fetchmany()
			if tups == []:
				break
			reslist.append(numpy.core.records.array(tups))
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
		raise Exception('Unrecognized driver')
	cur.close()
	con.commit()
	con.close()
	return res

def execute(query, db="wsdb", driver="psycopg2", user=None,
										password=None, host=None):
	if driver=='psycopg2':
		import pgdb
		con = pgdb.connect(database=db, user=user, password=password, host=host)
	elif driver=='sqlite3':
		import sqlite3
		con = sqlite3.connect(db)
	else:
		raise Exception('Unrecognized driver')
	cur=con.cursor()
	cur.execute(query)
	con.commit()
	con.close()

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
import time,pgdb


def fetchmany(curs, size=None, keep=False):
	"""Fetch the next set of rows of a query result.

	The number of rows to fetch per call is specified by the
	size parameter. If it is not given, the cursor's arraysize
	determines the number of rows to be fetched. If you set
	the keep parameter to true, this is kept as new arraysize.

	"""
	_cast = {'bool': numpy.bool,
		'int2': numpy.int16, 'int4': numpy.int32, 'serial': numpy.int32,
		'int8': numpy.int64, 'oid': numpy.int64, 'oid8': numpy.int64,
		'float4': numpy.float32, 'float8': numpy.float64,
		'numeric': numpy.float128, 'varchar': numpy.str}


	if size is None:
		size = curs.arraysize
	if keep:
		curs.arraysize = size
	result = curs._src.fetch(size)
	coltypes = [desc[1] for desc in curs.description]
	types=[_cast.get(a) for a in coltypes]
	return result,types #[list(row) for row in result]


def get(query, params=None, db="wsdb", driver="pgdb", user=None,
										password=None, host=None):
	'''This program executes the sql query and returns 
	the tuple of the numpy arrays.
	Example:
	a,b, c = sqlutil.get('select ra,dec,d25 from rc3')
	You can also use the parameters in your query:
	Example:
	a,b = squlil.get('select ra,dec from rc3 where name=?',"NGC 3166")
	'''
	if driver=='pgdb':
		import pgdb
		con = pgdb.connect(database=db, user=user, password=password, host=host)
	elif driver=='sqlite3':
		import sqlite3
		con = sqlite3.connect(db)
	else: 
		raise Exception("Unknown driver")
	cur = con.cursor()
	if params==None:
		res = cur.execute(query)
	else:
		res = cur.execute(query, params)

	if driver=='pgdb':
		tups,typelist = fetchmany(cur,-1, False)
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
	else:
		raise Exception('Unrecognized driver')
		

	cur.close()
	con.commit()
	con.close()
	nrows = len(tups)
#	print nrows
	if nrows == 0:
		return None
	for curtype in [numpy.float64,numpy.float32,numpy.int16,numpy.int32,numpy.int64]:
		if numpy.array([a ==curtype  for a in typelist]).all():
			res= numpy.asfarray(tups,typelist[0])
			return res.transpose()

	res = numpy.array(tups)
	nrows,ncols=res.shape
	result = []
	for  i, cur_type in zip(range(ncols),typelist):
		if cur_type == numpy.str:
			result.append(numpy.fromiter(res[:,i],dtype='|S100'))
			(result[-1])[result[-1]=='None']=''

		else:
			result.append(res[:,i].astype(cur_type))
	
	return result		
	#return [res[:,a].astype(b) for a,b in zip(range(ncols),types)]

def execute(query, db="wsdb", driver="pgdb", user=None,
										password=None, host=None):
	if driver=='pgdb':
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

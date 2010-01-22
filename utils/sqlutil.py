import types
import numpy
import time


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
	try:
		result = curs._src.fetch(size)
	except Error, err:
		raise pgdb.DatabaseError(str(err))
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
		con = sqlite3.connect(db, user=user, password=password, host=host)
	else: 
		raise Exception("Unknown driver")
	cur = con.cursor()
	if params==None:
		res = cur.execute(query)
	else:
		res = cur.execute(query, params)

	tups,types = fetchmany(cur,-1, False)

	cur.close()
	con.commit()
	con.close()
	nrows = len(tups)
#	print nrows
	if nrows == 0:
		return None
	for curtype in [numpy.float64,numpy.float32,numpy.int16,numpy.int32,numpy.int64]:
		if numpy.array([a ==curtype  for a in types]).all():
			res= numpy.asfarray(tups,types[0])
			return res.transpose()

	res = numpy.array(tups)
	nrows,ncols=res.shape
	result = []
	for  i, cur_type in zip(range(ncols),types):
		if cur_type == numpy.str:
			result.append(numpy.fromiter(res[:,i],dtype='|S100'))
			(result[-1])[result[-1]=='None']=''

		else:
			result.append(res[:,i].astype(cur_type))
	
	return result		
	#return [res[:,a].astype(b) for a,b in zip(range(ncols),types)]

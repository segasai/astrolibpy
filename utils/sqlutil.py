import types
import numpy

def get(query, params=None, db="wsdb", driver="pgdb"):
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
		con = pgdb.connect(database=db)
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
	tups = cur.fetchall()
	cur.close()
	con.commit()
	con.close()
	nrows = len(tups)
	if nrows == 0:
		return None
	line = tups[0]
	ncols = len(line)
	buf = []
	for i,curtype in enumerate(map(lambda s: type(s),line)):
		name = 'x%d'%i
		if curtype == types.StringType or curtype == types.UnicodeType:
			#TODO I'm now using the strings which are 100 char long 
			# probably it should be hardcoded that way.. 
			buf.append((name,'|S100'))
		else:
			buf.append((name,curtype))
	#doesn't really work because numpy.array wants the list of tuples, while
	#pgdb returns the list of list... 
	res = numpy.array(tups, dtype=numpy.dtype(buf))

	return [res['x%d'%i] for i in range(ncols)]

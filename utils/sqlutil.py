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
	res = numpy.array(tups)

	return res.transpose()

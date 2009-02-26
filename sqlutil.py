import pgdb
import sqlite3
import types
import numpy

def get(query, params=None, db="wsdb", driver="pgdb"):
	if driver=='pgdb':
		con = pgdb.connect(database=db)
	elif driver=='sqlite3':
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
	for curtype in map(lambda s: type(s),line):
		if curtype == types.StringType or curtype == types.UnicodeType:
			#TODO I'm now using the strings which are 100 char long 
			# probably it should be hardcoded that way.. 
			buf.append(numpy.zeros(nrows, dtype='|S100'))
		else:
			buf.append( numpy.zeros(nrows, dtype=curtype))

	for i in range(nrows):
		for j in range(ncols):
			buf[j][i]=tups[i][j]
	return tuple(buf)

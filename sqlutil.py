import pgdb
import types
import numpy

def get(query, db="wsdb"):
	con = pgdb.connect(database=db)
	cur = con.cursor()
	res = cur.execute(query)
	tups = cur.fetchall()
	nrows = len(tups)
	if nrows == 0:
		return None
	line = tups[0]
	ncols = len(line)
	buf = []
	for curtype in map(lambda s: type(s),line):
		if curtype == types.StringType:
			buf.append(range(nrows))
		else:
			buf.append( numpy.zeros(nrows, dtype=curtype))

	for i in range(nrows):
		for j in range(ncols):
			buf[j][i]=tups[i][j]
	return tuple(buf)

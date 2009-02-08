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
	buf = "("
	for i in range(ncols):
		curtype = type(line[i])
		if curtype == types.StringType:
			exec("a%d = range(nrows)" % i)
		else:
			exec("a%d = numpy.zeros(nrows, dtype=curtype)" % i)
		buf=buf+"a%d," % i
	buf = buf[:-1]+")"
	for i in range(nrows):
		for j in range(ncols):
			exec("a%d[i] = tups[i][j]" % j)
	#print buf
	#print a0
	exec("ret=" + buf)
	return ret
		
		


import scipy.io

def readcol(filename, separator=' ', format=None):
	if format==None:
		print 'x'
		res=scipy.io.read_array(filename, separator=separator)
		ncols = res.shape[1]
		nrows = res.shape[1]
		buf="("
		for i in range(ncols):
			exec("a%d = res[:,i]"%i)
			buf=buf+"a%d,"%i
		buf=buf[:-1]+")"
		exec("result=%s"%buf)
		return result
#	print 'z',format,format==None
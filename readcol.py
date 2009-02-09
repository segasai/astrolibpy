import scipy.io
import numpy 

def readcol(filename, delimiter=' ', format=None, skiprows=0):
	if format==None:
		print 'x'
		res=numpy.loadtxt(filename, delimiter=delimiter)
		ncols = res.shape[0]
		nrows = res.shape[1]
		stor=[]
		for i in range(ncols):
			stor.append(res[:,i])
		return tuple(stor)
	else:
		types=[]
		i=0
		formats=format.split(',')
		convs={}
		
		retnull = lambda s: numpy.float(s or 0)
		for i, a in enumerate(formats):
			if a=='I':
				curtype=numpy.int32
				convs[i]=retnull
			elif a=='F':
				curtype=numpy.float32
				convs[i]=retnull
			elif a=='D':
				curtype=numpy.float64
				convs[i]=retnull
			elif a=='S':
				curtype="S100"#numpy.str
			types.append(("a%d"%i,curtype))
		print types
		
		rec=numpy.loadtxt(file(filename),dtype=types, delimiter=delimiter, 
						skiprows=skiprows,converters=convs)
		print rec[0]
		ncols=len(rec[0])
		nrows=len(rec)

		buf="("
		stor=[]
		for a in formats:
			if a=='I':
				tmp=numpy.zeros(nrows,dtype=numpy.int32)
			elif a=='F':
				tmp=numpy.zeros(nrows,dtype=numpy.float32)
			elif a=='D':
				tmp=numpy.zeros(nrows,dtype=numpy.float64)
			elif a=='S':
				tmp=numpy.zeros(nrows,dtype=numpy.str)
			stor.append(tmp)

		for i in range(ncols):
			for j in range(nrows):
				stor[i][j]=rec[j][i]
		return tuple(stor)
#	print 'z',format,format==None

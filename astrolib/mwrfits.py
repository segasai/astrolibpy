import numpy, pyfits, types, itertools

def mwrfits(filename, arraylist, namelist=None, header=None):
	""" 
	Writes the list of numpy.arrays arraylist as a FITS table filename
	using namelist as list of names. 
	Arraylist can be dictionary with arrays as values and names as keys. 
	Also Arraylist can be numpy-record-array.
	Example:
	mwrfits('/tmp/xx.fits',[arr,arr1],['X','Y'])
	Or :
	mwrfits('test.fits',{'X':arr,'Y':arr1})
	Or:
	data = numpy.zeros((4,),dtype=[('run','i4'),('rerun','f8'),('zz','b')])
	mwfits('test1.fits',data)
	
	Keep in mind that when you used a dictionary, the order of columns in the
	fits file is not guaranteed
	"""
	tmplist=[]
	if isinstance(arraylist,numpy.ndarray):
		if arraylist.dtype.type is numpy.void:
			iter=itertools.izip(arraylist.dtype.names, itertools.imap (arraylist.__getitem__ , arraylist.dtype.names))
	else:
		if isinstance(arraylist,types.ListType) or isinstance(arraylist,types.TupleType):
			iter= zip(namelist, arraylist)
		elif isinstance(arraylist,types.DictType):
			iter= arraylist.iteritems()

	for name, arr in iter:
		if arr.dtype.type==numpy.int8:
			format='I'
		elif arr.dtype.type==numpy.int16:
			format='I'
		elif arr.dtype.type==numpy.int32:
			format='J'
		elif arr.dtype.type==numpy.int64:
			format='K'
		elif arr.dtype.type==numpy.float32:
			format='E'
		elif arr.dtype.type==numpy.float64:
			format='D'
		elif arr.dtype.type==numpy.string_:
			format='%dA'%arr.dtype.itemsize
		elif arr.dtype.type==numpy.bool_:
			format='I' # bool as int
		else:
			raise Exception("Oops unknown datatype %s"%arr.dtype)
		tmplist.append(pyfits.Column(name=name, array=arr, format=format))
	hdu = pyfits.new_table(tmplist)
	hdu.writeto(filename,clobber=True)
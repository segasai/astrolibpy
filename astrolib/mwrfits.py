import numpy, pyfits, types, itertools

def mwrfits(filename, arraylist, namelist=None):
	""" 
	Writes the list of numpy.arrays arraylist as a FITS table filename
	using namelist as list of names. 
	Arraylist can be dictionary with arrays as values and names as keys.
	Also Arraylist can be numpy-record-array
	"""
	tmplist=[]
	if isinstance(arraylist,numpy.ndarray):
		if arraylist.dtype.type is numpy.void:
			keys=arraylist.dtype.fields.iterkeys()
			iter=itertools.izip(keys, itertools.imap (arraylist.__getitem__ , keys))
	else:
		if isinstance(arraylist,types.ListType):
			iter= zip(namelist, arraylist)
		elif isinstance(arraylist,types.DictType):
			iter= arraylist.iteritems()

	for name, arr in iter:
		print name
		if arr.dtype==numpy.int16:
			format='I'
		elif arr.dtype==numpy.int32:
			format='J'
		elif arr.dtype==numpy.int64:
			format='K'
		elif arr.dtype==numpy.float32:
			format='E'
		elif arr.dtype==numpy.float64:
			format='D'
		else:
			raise Exception("Oops unknown datatype %s"%arr.dtype)
		tmplist.append(pyfits.Column(name=name, array=arr, format=format))
	hdu = pyfits.new_table(tmplist)
	hdu.writeto(filename,clobber=True)
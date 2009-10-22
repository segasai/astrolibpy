import numpy

def sphdist (ra1, dec1, ra2, dec2):
	"""measures the spherical distance in degrees
		The input has to be in degrees too
	"""
	dec1_r=numpy.deg2rad(dec1)
	dec2_r=numpy.deg2rad(dec2)
	return 2*numpy.rad2deg(numpy.arcsin(numpy.sqrt(
		(numpy.sin((dec1_r-dec2_r)/2))**2+ 
		numpy.cos(dec1_r)*numpy.cos(dec2_r)*
		(numpy.sin((numpy.deg2rad(ra1-ra2))/2))**2
		)))
	

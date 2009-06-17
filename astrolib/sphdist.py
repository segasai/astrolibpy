import numpy

def sphdist (ra1, dec1,ra2,dec2):
	"""measures the spherical distance in degrees"""
	ra1r=numpy.deg2rad(ra1)
	ra2r=numpy.deg2rad(ra2)
	dec1r=numpy.deg2rad(dec1)
	dec2r=numpy.deg2rad(dec2)
	return 2*numpy.rad2deg(numpy.arcsin(numpy.sqrt(
		(numpy.sin((dec1r-dec2r)/2))**2+ 
		numpy.cos(dec1r)*numpy.cos(dec2r)*(numpy.sin((ra1r-ra2r)/2))**2
		)))
	

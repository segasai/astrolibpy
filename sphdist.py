import numpy

def sphdist (ra1, dec1,ra2,dec2):
	"""measures the spherical distance in degrees"""
	fac=180./numpy.pi
	ra1r=ra1/fac
	ra2r=ra2/fac
	dec1r=dec1/fac
	dec2r=dec2/fac
	return 2*fac*numpy.arcsin(numpy.sqrt(
		(numpy.sin((dec1r-dec2r)/2))**2+ 
		numpy.cos(dec1r)*numpy.cos(dec2r)*(numpy.sin((ra1r-ra2r)/2))**2
		))
	

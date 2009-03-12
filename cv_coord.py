import numpy
def cv_coord(a,b,c,fr=None,to=None):
	radec = 180./numpy.pi
	if fr=='sph':
		cosa = numpy.cos(a/radec)
		sina = numpy.sin(a/radec)
		cosb = numpy.cos(b/radec)
		sinb = numpy.sin(b/radec)
		x=c*cosa*cosb
		y=c*sina*cosb
		z=c*sinb
	if fr=='rect':
		x=a
		y=b
		z=c
	if to=='rect':
		return (x,y,z)
	if to=='sph':
		ra = numpy.arctan2(y,x)*radec
		dec = numpy.arctan2(z,numpy.sqrt(x**2+y**2))*radec
		rad = numpy.sqrt(x**2+y**2+z**2)
		return (ra,dec,rad)
import numpy
def cv_coord(a,b,c,fr=None,to=None,degr=False):
	if degr:
		radec = 180./numpy.pi
	else:
		radec= 1
	if fr=='sph':
		cosa = numpy.cos(a/radec)
		sina = numpy.sin(a/radec)
		cosb = numpy.cos(b/radec)
		sinb = numpy.sin(b/radec)
		x=c*cosa*cosb
		y=c*sina*cosb
		z=c*sinb
	elif fr=='rect':
		x=a
		y=b
		z=c
	elif fr is None:
		raise Exception('Please specify the input coordinate system')
	else:
		raise Exception('Unknown input coordinate system')
	if to=='rect':
		return (x,y,z)
	elif to=='sph':
		ra = numpy.arctan2(y,x)*radec
		dec = numpy.arctan2(z,numpy.sqrt(x**2+y**2))*radec
		rad = numpy.sqrt(x**2+y**2+z**2)
		return (ra,dec,rad)
	elif to is None:
		raise Exception('Please specify the output coordinate system')
	else:
		raise Exception('Unknown output coordinate system')

		
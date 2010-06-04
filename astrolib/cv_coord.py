import numpy

def cv_coord(a,b,c,fr=None,to=None,degr=False):
	if degr:
		degrad = numpy.deg2rad
		raddeg = numpy.rad2deg
	else:
		degrad = lambda x: x
		raddeg = lambda x: x
	if fr=='sph':
		cosa = numpy.cos(degrad(a))
		sina = numpy.sin(degrad(a))
		cosb = numpy.cos(degrad(b))
		sinb = numpy.sin(degrad(b))
		x=c*cosa*cosb
		y=c*sina*cosb
		z=c*sinb
	elif fr=='rect':
		x=a
		y=b
		z=c
	elif fr is None:
		raise Exception('You must specify the input coordinate system')
	else:
		raise Exception('Unknown input coordinate system')
	if to=='rect':
		return (x,y,z)
	elif to=='sph':
		ra = raddeg(numpy.arctan2(y,x))
		dec = raddeg(numpy.arctan2(z,numpy.sqrt(x**2+y**2)))
		rad = numpy.sqrt(x**2+y**2+z**2)
		return (ra,dec,rad)
	elif to is None:
		raise Exception('You must specify the output coordinate system')
	else:
		raise Exception('Unknown output coordinate system')

		
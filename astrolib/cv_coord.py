from numpy import sin,cos,deg2rad,rad2deg,arctan2,sqrt

def cv_coord(a,b,c,fr=None,to=None,degr=False):
	if degr:
		degrad = deg2rad
		raddeg = rad2deg
	else:
		degrad = lambda x: x
		raddeg = lambda x: x
	if fr=='sph':
		x=c*cos(degrad(a))*cos(degrad(b))
		y=c*sin(degrad(a))*cos(degrad(b))
		z=c*sin(degrad(b))
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
		ra = raddeg(arctan2(y,x))
		dec = raddeg(arctan2(z,sqrt(x**2+y**2)))
		rad = sqrt(x**2+y**2+z**2)
		return (ra,dec,rad)
	elif to is None:
		raise Exception('You must specify the output coordinate system')
	else:
		raise Exception('Unknown output coordinate system')

		
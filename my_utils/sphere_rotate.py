import numpy
import numexpr
from cv_coord import cv_coord

def torect(ra,dec):
	x=numexpr.evaluate('cos(ra/57.295779513082323)*cos(dec/57.295779513082323)')
	y=numexpr.evaluate('sin(ra/57.295779513082323)*cos(dec/57.295779513082323)')
	z=numexpr.evaluate('sin(dec/57.295779513082323)')
	return x,y,z

def fromrect(x,y,z):
	ra=numexpr.evaluate('arctan2(y,x)*57.295779513082323')
	dec=numexpr.evaluate('57.295779513082323*arctan2(z,sqrt(x**2+y**2))')
	return ra,dec

def sphere_rotate(ra, dec, rapol, decpol, ra0, revert=False):
	""" rotate ra,dec to a new spherical coordinate system where the pole is 
		at rapol,decpol and the zeropoint is at ra=ra0 
		revert flag allows to reverse the transformation
	"""

	x,y,z=torect(ra,dec)

	tmppol=cv_coord(rapol,decpol,1,degr=True,fr='sph',to='rect') #pole axis
	tmpvec1=cv_coord(ra0,0,1,degr=True,fr='sph',to='rect') #x axis
	tmpvec1=numpy.array(tmpvec1)

	tmpvec1[2]=(-tmppol[0]*tmpvec1[0]-tmppol[1]*tmpvec1[1])/tmppol[2]
	tmpvec1/=numpy.sqrt((tmpvec1**2).sum())
	tmpvec2=numpy.cross(tmppol,tmpvec1) # y axis 

	if not revert:
		Axx,Axy,Axz=tmpvec1
		Ayx,Ayy,Ayz=tmpvec2
		Azx,Azy,Azz=tmppol
	else:
		Axx,Ayx,Azx=tmpvec1
		Axy,Ayy,Azy=tmpvec2
		Axz,Ayz,Azz=tmppol
	xnew = numexpr.evaluate('x*Axx+y*Axy+z*Axz')
	ynew = numexpr.evaluate('x*Ayx+y*Ayy+z*Ayz')
	znew = numexpr.evaluate('x*Azx+y*Azy+z*Azz')
	
	del x,y,z
	tmp = fromrect(xnew,ynew,znew)
	return (tmp[0],tmp[1])

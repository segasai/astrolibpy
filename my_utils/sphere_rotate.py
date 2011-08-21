import numpy
from cv_coord import cv_coord


def sphere_rotate(ra, dec, rapol, decpol, ra0, revert=False):
	""" rotate ra,dec to a new spherical coordinate system where the pole is 
		at rapol,decpol and the zeropoint is at ra=ra0 
		revert flag allows to reverse the transformation
	"""
	x,y,z=cv_coord(ra,dec,dec*0+1,fr='sph',to='rect',degr=True)

	tmppol=cv_coord(rapol,decpol,1,degr=True,fr='sph',to='rect') #pole axis
	tmpvec1=cv_coord(ra0,0,1,degr=True,fr='sph',to='rect') #x axis
	tmpvec1=numpy.array(tmpvec1)

	tmpvec1[2]=(-tmppol[0]*tmpvec1[0]-tmppol[1]*tmpvec1[1])/tmppol[2]
	tmpvec1/=numpy.sqrt((tmpvec1**2).sum())
	tmpvec2=numpy.cross(tmppol,tmpvec1) # y axis 

	if not revert:
		xnew=x*tmpvec1[0]+y*tmpvec1[1]+z*tmpvec1[2]
		ynew=x*tmpvec2[0]+y*tmpvec2[1]+z*tmpvec2[2]
		znew=x*tmppol[0] +y*tmppol[1] +z*tmppol[2]
	else:
		xnew=x*tmpvec1[0]+y*tmpvec2[0]+z*tmppol[0]
		ynew=x*tmpvec1[1]+y*tmpvec2[1]+z*tmppol[1]
		znew=x*tmpvec1[2]+y*tmpvec2[2]+z*tmppol[2]
	del x,y,z
	tmp=cv_coord(xnew,ynew,znew,to='sph',fr='rect',degr=True)
	return (tmp[0],tmp[1])

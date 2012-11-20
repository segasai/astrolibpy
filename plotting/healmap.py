# coding: utf-8
import numpy,healpy,matplotlib.collections,matplotlib.patches
import sqlutil,quick_hist
import matplotlib.pyplot as plt 
from idlplot import plot,oplot

import healpy
import scipy as sc

# Joe Bovy's code of his extensions of the healpy module
# Copyright Joe Bovy
class bovy_healpy:
	@staticmethod
	def pix2vert(nside,ipix,nest=False):
		"""
		NAME:
		   pix2vert
		PURPOSE:
		   calculate the locations of the vertices (theta,phi)
		   of a given HEALPix pixel
		INPUT:
		   nside - HEALPix resolution parameter
		   ipix - pixel number
		   nest - if True, use NESTED scheme (default: RING)
		OUTPUT:
		   numpy.array([4,2]) theta,phi [rad] NWSE
		HISTORY:
		   2010-01-21 - Written - Bovy (NYU)
		"""
		(centerTheta,centerPhi)= healpy.pix2ang(nside,ipix,nest=nest)
		#Are we in the polar regime or in the equatorial regime?
		z= sc.cos(centerTheta)
		if z > -2./3. and z < 2./3.:
			return bovy_healpy._ang2vert_eq(nside,centerTheta,centerPhi,z)
		else:
			return bovy_healpy._ang2vert(nside,centerTheta,centerPhi,z)
	@staticmethod		
	def _ang2vert_eq(nside,theta,phi,z):
		a= 4./3./nside
		b= 8./3./sc.pi
		deltaZ= a
		deltaPhi= a/b
		out= []
		out.append([sc.arccos(z+deltaZ/2),phi])
		out.append([theta,phi-deltaPhi/2.])
		out.append([sc.arccos(z-deltaZ/2),phi])
		out.append([theta,phi+deltaPhi/2.])
		return sc.array(out)
	@staticmethod
	def _ang2vert(nside,theta,phi,z):
		(xsCenter,ysCenter)= bovy_healpy._ang2xsys(z,phi)
		delta= sc.pi/4./nside
		out= []
		out.append(bovy_healpy._xsys2ang(xsCenter,ysCenter+delta))
		out.append(bovy_healpy._xsys2ang(xsCenter-delta,ysCenter))
		out.append(bovy_healpy._xsys2ang(xsCenter,ysCenter-delta))
		out.append(bovy_healpy._xsys2ang(xsCenter+delta,ysCenter))
		return sc.array(out)
	@staticmethod
	def _xsys2ang(xs,ys):
		if sc.fabs(ys) <= sc.pi/4.:
			return [sc.arccos(8./3./sc.pi*ys),xs]
		else:
			xt= (xs % (sc.pi/2.))
			fabsys= sc.fabs(ys)
			theta= sc.arccos((1.-1./3.*(2.-4.*fabsys/sc.pi)**2.)*ys/fabsys)
			if fabsys == sc.pi/2.:
				phi= xs-fabsys+sc.pi/4.
			else:
				phi= xs-(fabsys-sc.pi/4.)/(fabsys-sc.pi/2.)*(xt-sc.pi/4.) 
			return [theta % (sc.pi+0.0000000001),phi % (2.*sc.pi)] #Hack
	@staticmethod
	def _ang2xsys(z,phi):
		if sc.fabs(z) <= 2./3.:
			return [phi,3.*sc.pi/8.*z]
		else:
			phit= (phi % (sc.pi/2.))
			sigz= bovy_healpy._sigma(z)
			return [phi-(sc.fabs(sigz)-1.)*(phit-sc.pi/4.),sc.pi/4.*sigz]

	@staticmethod
	def _sigma(z):
		if z < 0.:
			return -(2.-sc.sqrt((3.*(1+z))))

		else:
			return 2.-sc.sqrt((3.*(1-z)))
					
def healmap(ras, decs, ramin=0, ramax=360, decmin=-90, decmax=90, nside=64,
			xflip=False, vmax=None, vmin=None, cmap=plt.cm.gray_r,
			linewidth=1, weights=None, wrap_angle=None, skip_empty=True,
			weight_norm=False, rasterized=True):
	"""
		Make the 2D histogram of the datapoints on the sky using the Healpix
		pixels 
	"""
	
	plt.ioff()

	nel = 12 * nside**2
	ipix = healpy.ang2pix(nside, numpy.deg2rad(decs)+numpy.pi*0.5, numpy.deg2rad(ras))
	hh = quick_hist.quick_hist((ipix,), nbins=[nel],
		range=[[0, nel]],weights=weights)
	hhsum = quick_hist.quick_hist((ipix,), nbins=[nel],
		range=[[0, nel]])

	pixarea = (4*numpy.pi*(180/numpy.pi)**2)/nel # in sq. deg
	hh = hh / pixarea
	hhsum = hhsum / pixarea

	xloc = numpy.arange(nel)
	verts = [bovy_healpy.pix2vert(nside, xx) for xx in xloc]

	collist = []
	fac = 180/numpy.pi
	if wrap_angle is None:
		wrap_angle = 2 * numpy.pi
	else:
		wrap_angle = numpy.deg2rad(wrap_angle)
	for ii, v in enumerate(verts):
		if skip_empty and hhsum[ii] == 0:
			continue
		
		b,a = v.T
		if (a.max()-a.min())>(numpy.pi):
			a = ((a + (numpy.pi - a[0])) % (2 * numpy.pi)) - (numpy.pi-a[0])
		if a.mean() < 0:
			a = (a + 2 * numpy.pi)
		if a.mean() > wrap_angle:
			a -= 2 * numpy.pi
		xys = numpy.array([a * fac, b * fac - 90]).T

		curpoly = matplotlib.patches.Polygon(xys, closed=True)
		collist.append(curpoly)

	raranges = [ramin,ramax]
	if xflip:
		raranges = raranges[::-1]
	plot([-1], xr=raranges, yr=[decmin,decmax])
	ax = plt.gca()
	tmppol = matplotlib.patches.Polygon(numpy.array(([0, 360, 360, 0],
										[-90, -90, 90, 90])).T,closed=True)
	tmpcoll = matplotlib.collections.PatchCollection([tmppol], linewidths=0.01,
		edgecolors='black', facecolors='black')
	tmpcoll.set_array(numpy.array([0]))


	coll = matplotlib.collections.PatchCollection(collist,
		linewidths=linewidth)

	if vmin is None:
		vmin = 0 
	if vmax is None:
		vmax = hh.max()

	if weight_norm:
		hh = hh * 1. / (hhsum + 1 * (hhsum == 0))
	
	if skip_empty:
		hh1 = hh[hhsum > 0]
	else:
		hh1 = hh

	coll.set_array(hh1)
	coll.set_cmap(cmap)
	coll.set_clim(vmin, vmax)
	coll.set_edgecolors(coll.cmap(coll.norm(hh1)))
	coll.set_rasterized(rasterized)
	tmpcoll.set_clim(vmin, vmax)

	ax.add_collection(coll)
	return coll
	
	
if __name__=='__main__':
	ras,decs=sqlutil.get('select radeg,dedeg from sdss_phot_qso.main',
		host='cappc118')
	doit(ras,decs, ramin=0,ramax=360,decmin=-90, decmax=90, nside=64)

	plt.savefig('xx.png',dpi=200)

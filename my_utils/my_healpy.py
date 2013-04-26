import healpy,numpy as np

def radec2pix(nside,ra,dec,nest=True):
	_,__=np.pi/2-np.deg2rad(dec),np.deg2rad(ra)
	return healpy.ang2pix(nside,_,__,nest=nest)

def pix2radec(nside,pix,nest=True):
	__,_ = healpy.pix2ang(nside,pix,nest=nest)
	ra,dec=np.rad2deg(_),np.rad2deg(np.pi/2-__)
	return ra,dec

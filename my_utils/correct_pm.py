import gal_uvw, cv_coord, euler
import numpy as np

def get_uvw_sun(vlsr=220):
	
	usun, vsun, wsun = -8.5, 13.38 + vlsr, 6.49 # the signs are in the coord system of gal_uvw
											 # e.g. U points toward anticenter
	return usun,vsun,wsun

def correct_pm(ra, dec, pmra, pmdec, dist, vlsr=220):
	"""Corrects the proper motion for the speed of the Sun
	Arguments:
		ra - RA in deg
		dec -- Declination in deg
		pmra -- pm in RA in mas/yr
		pmdec -- pm in declination in mas/yr
		dist -- distance in kpc
	Returns:
		(pmra,pmdec) the tuple with the proper motions corrected for the Sun's motion
	"""
	
	one = ra * 0 + 1
	zero = ra * 0
	usun,vsun,wsun = get_uvw_sun(vlsr=vlsr)
	dist_pc = dist * 1000. 
	ur, vr, wr = gal_uvw.gal_uvw(distance=dist_pc, ra=ra, dec=dec,
		pmra=one, pmdec=zero, vrad=zero)
	ud, vd, wd = gal_uvw.gal_uvw(distance=dist_pc, ra=ra, dec=dec,
		pmra=zero, pmdec=one, vrad=zero)



	d_pmra =  -(ur * usun + vr * vsun + wr * wsun) / (ur**2 + vr**2 + wr**2)
	d_pmdec = -(ud * usun + vd * vsun + wd * wsun) / (ud**2 + vd**2 + wd**2)
	# d_pmra d_pmdec -- these should be the pm's of the non-moving object as seen from the moving Sun
	return (pmra-d_pmra,pmdec-d_pmdec)


def correct_vel(ra, dec, vel, vlsr=220):
	"""Corrects the proper motion for the speed of the Sun
	Arguments:
		ra - RA in deg
		dec -- Declination in deg
		pmra -- pm in RA in mas/yr
		pmdec -- pm in declination in mas/yr
		dist -- distance in kpc
	Returns:
		(pmra,pmdec) the tuple with the proper motions corrected for the Sun's motion
	"""
	
	
	l,b = euler.euler(ra, dec)
	l = np.deg2rad(l)
	b = np.deg2rad(b)
	usun,vsun,wsun = get_uvw_sun(vlsr=vlsr)

	delta = -usun*np.cos(l)*np.cos(b) + vsun * np.sin(l) * np.cos(b) + wsun * np.sin(b)
	# projection of the sun's velocity to the lign of sight vector
	# notice the first minus -- it is because the usun is in the coord system where 
	# X points towards anticenter

	#return the corrected velocity
	return vel + delta

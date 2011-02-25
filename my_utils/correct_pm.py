import gal_uvw, cv_coord, euler


def correct_pm(ra, dec, pmra, pmdec, dist):
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

	dist_pc = dist * 1000. 
	ur, vr, wr = gal_uvw.gal_uvw(distance=dist_pc, ra=ra, dec=dec,
		pmra=one, pmdec=zero, vrad=zero)
	ud, vd, wd = gal_uvw.gal_uvw(distance=dist_pc, ra=ra, dec=dec,
		pmra=zero, pmdec=one, vrad=zero)

	vlsr = 220


	usun, vsun, wsun = -8.5, 13.38 + vlsr, 6.49 # the signs are in the coord system of gal_uvw
												 # e.g. U points toward anticenter

	d_pmra =  -(ur * usun + vr * vsun + wr * wsun) / (ur**2 + vr**2 + wr**2)
	d_pmdec = -(ud * usun + vd * vsun + wd * wsun) / (ud**2 + vd**2 + wd**2)
	# d_pmra d_pmdec -- these should be the pm's of the non-moving object as seen from the moving Sun
	return (pmra-d_pmra,pmdec-d_pmdec)

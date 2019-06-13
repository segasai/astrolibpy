import gal_uvw, cv_coord, euler
import numpy as np
import astropy.coordinates as acoo
import astropy.units as auni

def get_uvw_sun(vlsr):
    
    usun, vsun, wsun = -11.1, 12.24 + vlsr, 7.25 # the signs are in the coord system of gal_uvw
	# this is from schonrich, binney        
	# e.g. U points toward anticenter
    return usun,vsun,wsun

vlsr0 = 232.8 # from mcmillan 2017

def correct_pm(ra, dec, pmra, pmdec, dist, vlsr=vlsr0):
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
    usun,vsun,wsun = get_uvw_sun(vlsr)
    dist_pc = dist * 1000. 
    ur, vr, wr = gal_uvw.gal_uvw(distance=dist_pc, ra=ra, dec=dec,
        pmra=one, pmdec=zero, vrad=zero)
    ud, vd, wd = gal_uvw.gal_uvw(distance=dist_pc, ra=ra, dec=dec,
        pmra=zero, pmdec=one, vrad=zero)



    d_pmra =  -(ur * usun + vr * vsun + wr * wsun) / (ur**2 + vr**2 + wr**2)
    d_pmdec = -(ud * usun + vd * vsun + wd * wsun) / (ud**2 + vd**2 + wd**2)
    # d_pmra d_pmdec -- these should be the pm's of the non-moving object as seen from the moving Sun
    return (pmra-d_pmra,pmdec-d_pmdec)


def correct_vel(ra, dec, vel, vlsr=vlsr0):
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
    
    C=acoo.ICRS(ra=ra*auni.deg,dec=dec*auni.deg,
                    radial_velocity=vel*auni.km/auni.s,
                    distance=np.ones_like(vel)*auni.kpc,
                    pm_ra_cosdec=np.zeros_like(vel)*auni.mas/auni.year,
                    pm_dec=np.zeros_like(vel)*auni.mas/auni.year)
    #frame = acoo.Galactocentric (galcen_vsun = np.array([ 11.1, vlsr+12.24, 7.25])*auni.km/auni.s)
    kw = dict(galcen_v_sun = acoo.CartesianDifferential(np.array([ 11.1, vlsr+12.24, 7.25])*auni.km/auni.s))                                 
    frame = acoo.Galactocentric (**kw)
    Cg = C.transform_to(frame)
    Cg1 = acoo.Galactocentric(x=Cg.x, y=Cg.y, z=Cg.z,
                              v_x=Cg.v_x*0,
                              v_y=Cg.v_y*0,
                              v_z=Cg.v_z*0, **kw)
    C1=Cg1.transform_to(acoo.ICRS)
    return np.asarray(((C.radial_velocity-C1.radial_velocity)/(auni.km/auni.s)).decompose())

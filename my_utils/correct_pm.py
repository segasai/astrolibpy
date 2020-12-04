import numpy as np
import astropy.coordinates as acoo
import astropy.units as auni

vlsr0 = 232.8  # from mcmillan 2017


def correct_pm(ra, dec, pmra, pmdec, dist, vlsr=vlsr0, split=None):
    if split is None:
        return correct_pm0(ra, dec, pmra, pmdec, dist, vlsr=vlsr0)
    else:
        N = len(ra)
        n1 = N // split

        ra1 = np.array_split(ra, n1)
        dec1 = np.array_split(dec, n1)
        pmra1 = np.array_split(pmra, n1)
        pmdec1 = np.array_split(pmdec, n1)
        if hasattr(dist, '__len__'):
            assert (len(dist) == N)
        else:
            dist = np.zeros(N) + dist
        dist1 = np.array_split(dist, n1)
        ret = []
        for curra, curdec, curpmra, curpmdec, curdist in zip(
                ra1, dec1, pmra1, pmdec1, dist1):
            ret.append(
                correct_pm0(curra,
                            curdec,
                            curpmra,
                            curpmdec,
                            curdist,
                            vlsr=vlsr0))
        retpm1 = np.concatenate([_[0] for _ in ret])
        retpm2 = np.concatenate([_[1] for _ in ret])
        return retpm1, retpm2


def correct_pm0(ra, dec, pmra, pmdec, dist, vlsr=vlsr0):
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
    C = acoo.ICRS(ra=ra * auni.deg,
                  dec=dec * auni.deg,
                  radial_velocity=0 * auni.km / auni.s,
                  distance=dist * auni.kpc,
                  pm_ra_cosdec=pmra * auni.mas / auni.year,
                  pm_dec=pmdec * auni.mas / auni.year)
    kw = dict(galcen_v_sun=acoo.CartesianDifferential(
        np.array([11.1, vlsr + 12.24, 7.25]) * auni.km / auni.s))
    frame = acoo.Galactocentric(**kw)
    Cg = C.transform_to(frame)
    Cg1 = acoo.Galactocentric(x=Cg.x,
                              y=Cg.y,
                              z=Cg.z,
                              v_x=Cg.v_x * 0,
                              v_y=Cg.v_y * 0,
                              v_z=Cg.v_z * 0,
                              **kw)
    C1 = Cg1.transform_to(acoo.ICRS())
    return ((C.pm_ra_cosdec - C1.pm_ra_cosdec).to_value(auni.mas / auni.year),
            (C.pm_dec - C1.pm_dec).to_value(auni.mas / auni.year))


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

    C = acoo.ICRS(ra=ra * auni.deg,
                  dec=dec * auni.deg,
                  radial_velocity=vel * auni.km / auni.s,
                  distance=np.ones_like(vel) * auni.kpc,
                  pm_ra_cosdec=np.zeros_like(vel) * auni.mas / auni.year,
                  pm_dec=np.zeros_like(vel) * auni.mas / auni.year)
    kw = dict(galcen_v_sun=acoo.CartesianDifferential(
        np.array([11.1, vlsr + 12.24, 7.25]) * auni.km / auni.s))
    frame = acoo.Galactocentric(**kw)
    Cg = C.transform_to(frame)
    Cg1 = acoo.Galactocentric(x=Cg.x,
                              y=Cg.y,
                              z=Cg.z,
                              v_x=Cg.v_x * 0,
                              v_y=Cg.v_y * 0,
                              v_z=Cg.v_z * 0,
                              **kw)
    C1 = Cg1.transform_to(acoo.ICRS())
    return np.asarray(((C.radial_velocity - C1.radial_velocity) /
                       (auni.km / auni.s)).decompose())

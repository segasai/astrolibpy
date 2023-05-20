import numpy as np
import astropy.coordinates as acoo
import astropy.units as auni

# Use v4.0 defaults
GCPARAMS = acoo.galactocentric_frame_defaults.get_from_registry(
    "v4.0")['parameters']

kms = auni.km / auni.s
masyr = auni.mas / auni.year


def correct_pm(ra, dec, pmra, pmdec, dist, split=None, vlsr=None):
    """Corrects the proper motion for the speed of the Sun
    Arguments:
        ra - RA in deg
        dec -- Declination in deg
        pmra -- pm in RA in mas/yr (with cosine term)
        pmdec -- pm in declination in mas/yr
        dist -- distance in kpc
        split (optional) -- integer number to split calculations in blocks
               to avoid OOM
    Returns:
        (pmra,pmdec) the tuple with the proper motions corrected for the
        Sun's motion
    """
    if vlsr is not None:
        print('WARNING vlsr is ignored')
    if split is None:
        return _correct_pm0(ra, dec, pmra, pmdec, dist)
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
            ret.append(_correct_pm0(curra, curdec, curpmra, curpmdec, curdist))
        retpm1 = np.concatenate([_[0] for _ in ret])
        retpm2 = np.concatenate([_[1] for _ in ret])
        return retpm1, retpm2


def _correct_pm0(ra, dec, pmra, pmdec, dist):
    """Corrects the proper motion for the speed of the Sun
    Arguments:
        ra - RA in deg
        dec -- Declination in deg
        pmra -- pm in RA in mas/yr (with a cosine term)
        pmdec -- pm in declination in mas/yr
        dist -- distance in kpc
    Returns:
        (pmra,pmdec) the tuple with the proper motions corrected for the
        Sun's motion
    """
    C = acoo.ICRS(ra=ra * auni.deg,
                  dec=dec * auni.deg,
                  radial_velocity=0 * kms,
                  distance=dist * auni.kpc,
                  pm_ra_cosdec=pmra * masyr,
                  pm_dec=pmdec * masyr)
    frame = acoo.Galactocentric(**GCPARAMS)
    Cg = C.transform_to(frame)
    Cg1 = acoo.Galactocentric(x=Cg.x,
                              y=Cg.y,
                              z=Cg.z,
                              v_x=Cg.v_x * 0,
                              v_y=Cg.v_y * 0,
                              v_z=Cg.v_z * 0,
                              **GCPARAMS)
    C1 = Cg1.transform_to(acoo.ICRS())
    return ((C.pm_ra_cosdec - C1.pm_ra_cosdec).to_value(masyr),
            (C.pm_dec - C1.pm_dec).to_value(masyr))


def correct_vel(ra, dec, vel, vlsr=None, split=None):
    """Corrects the proper motion for the speed of the Sun
    Arguments:
        ra - RA in deg
        dec -- Declination in deg
        vel -- heliocentric rv
    Returns:
        radial velocity corrected for solar reflex motion
    """
    if vlsr is not None:
        print('WARNING vlsr is ignored')
    if split is None:
        return _correct_vel0(ra, dec, vel)
    else:
        N = len(ra)
        n1 = N // split

        ra1 = np.array_split(ra, n1)
        dec1 = np.array_split(dec, n1)
        vel1 = np.array_split(vel, n1)
        ret = []
        for curra, curdec, curvel in zip(ra1, dec1, vel1):
            ret.append(_correct_vel0(curra, curdec, curvel))
        ret = np.concatenate(ret)
        return ret


def _correct_vel0(ra, dec, vel, vlsr=None):
    """Corrects the proper motion for the speed of the Sun
    Arguments:
        ra - RA in deg
        dec -- Declination in deg
        pmra -- pm in RA in mas/yr
        pmdec -- pm in declination in mas/yr
        dist -- distance in kpc
    Returns:
        (pmra,pmdec) the tuple with the proper motions corrected for the
        Sun's motion
    """
    if vlsr is not None:
        print('WARNING vlsr is ignored')

    C = acoo.ICRS(ra=ra * auni.deg,
                  dec=dec * auni.deg,
                  radial_velocity=vel * kms,
                  distance=np.ones_like(vel) * auni.kpc,
                  pm_ra_cosdec=np.zeros_like(vel) * masyr,
                  pm_dec=np.zeros_like(vel) * masyr)
    frame = acoo.Galactocentric(**GCPARAMS)
    Cg = C.transform_to(frame)
    Cg1 = acoo.Galactocentric(x=Cg.x,
                              y=Cg.y,
                              z=Cg.z,
                              v_x=Cg.v_x * 0,
                              v_y=Cg.v_y * 0,
                              v_z=Cg.v_z * 0,
                              **GCPARAMS)
    C1 = Cg1.transform_to(acoo.ICRS())
    return np.asarray((C.radial_velocity - C1.radial_velocity).to_value(kms))

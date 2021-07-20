import sphere_rotate
import numpy as np


def rotate_pm(ra, dec, pmra, pmdec, rapol, decpol, ra0, revert=False):
    """ 
    Rotate the proper motion to the sphere_rotate coord system 
    using finite differences. 
    I assume all angles are indeg and proper motions are in mas/yr
    Arguments:
    ra: 
    dec:
    pmra:
    pmdec:
    rapol: float
        RA of the pole
    decpol float
        Dec of the pole
    ra0: float
        RA of the (0,0) point of the new coordinate system 
    Returns:
    pmphi1, pmphi2 in the new coordinate system
    """

    fi1, fi2 = sphere_rotate.sphere_rotate(ra,
                                           dec,
                                           rapol,
                                           decpol,
                                           ra0,
                                           revert=revert)
    dt = 1  # year
    mult = 3600e3  # conversion from  degrees to mas
    ra1 = ra + pmra * dt / mult / np.cos(np.deg2rad(dec))
    dec1 = dec + pmdec * dt / mult
    fi1_1, fi2_1 = sphere_rotate.sphere_rotate(ra1,
                                               dec1,
                                               rapol,
                                               decpol,
                                               ra0,
                                               revert=revert)
    del ra1, dec1
    dpm1 = (fi1_1 - fi1) * np.cos(np.deg2rad(fi2)) * mult / dt
    dpm2 = (fi2_1 - fi2) * mult / dt
    return dpm1, dpm2

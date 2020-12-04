import sphere_rotate
import numpy as np


def rotate_pm(ra, dec, pmra, pmdec, rapol, decpol, ra0):
    # rotate the proper motion to the sphere_rotate coord sys
    fi1, fi2 = sphere_rotate.sphere_rotate(ra, dec, rapol, decpol, ra0)
    dt = 1  # year
    mult = 3600e3
    ra1 = ra + pmra * dt / mult / np.cos(np.deg2rad(dec))
    dec1 = dec + pmdec * dt / mult
    fi1_1, fi2_1 = sphere_rotate.sphere_rotate(ra1, dec1, rapol, decpol, ra0)
    del ra1, dec1
    dpm1 = (fi1_1 - fi1) * np.cos(np.deg2rad(fi2)) * mult / dt
    dpm2 = (fi2_1 - fi2) * mult / dt
    return dpm1, dpm2

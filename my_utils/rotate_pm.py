import sphere_rotate
import numpy as np


def cosd(x):
    return np.cos(np.deg2rad(x))


def sind(x):
    return np.sin(np.deg2rad(x))


def rotate_pm(ra, dec, pmra, pmdec, rapol, decpol, ra0, revert=False):
    """
    Rotate the proper motion to the sphere_rotate coord system 
    that is specified by the pole direction and right ascencion of the (0,0) pt
    I assume all angles are in degrees and proper motions are in mas/yr
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
    revert: bool
        if true do the inverse transformation. I.e. convert 
        phi1,phi2,pmphii1,pmphi2
        to pmra,pmdec
    Returns:
    pmphi1, pmphi2 in the new coordinate system
    """

    ra, dec, pmra, pmdec = [np.atleast_1d(_) for _ in [ra, dec, pmra, pmdec]]
    M = sphere_rotate.rotation_matrix(rapol, decpol, ra0)
    if revert:
        M = M.T
    # unit vectors
    e_mura = np.array([-sind(ra), cosd(ra), ra * 0])
    e_mudec = np.array(
        [-sind(dec) * cosd(ra), -sind(dec) * sind(ra),
         cosd(dec)])
    # velocity vector in arbitrary units
    V = pmra * e_mura + pmdec * e_mudec
    del e_mura, e_mudec
    # apply rotation to velocity
    V1 = M @ V
    del V
    X = np.array([cosd(ra) * cosd(dec), sind(ra) * cosd(dec), sind(dec)])
    # apply rotation to position
    X1 = M @ X
    del X
    # rotated coordinates in radians
    lon = np.arctan2(X1[1, :], X1[0, :])
    lat = np.arctan2(X1[2, :], np.sqrt(X1[0, :]**2 + X1[1, :]**2))
    del X1
    # unit vectors in rotated coordinates
    e_mura = np.array([-np.sin(lon), np.cos(lon), lon * 0])
    e_mudec = np.array(
        [-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon),
         np.cos(lat)])
    del lon, lat
    return np.sum(e_mura * V1, axis=0), np.sum(e_mudec * V1, axis=0)

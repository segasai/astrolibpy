import numpy as np
import numexpr


def torect(ra, dec):
    x = numexpr.evaluate(
        'cos(ra/57.295779513082323)*cos(dec/57.295779513082323)')
    y = numexpr.evaluate(
        'sin(ra/57.295779513082323)*cos(dec/57.295779513082323)')
    z = numexpr.evaluate('sin(dec/57.295779513082323)')
    return x, y, z


def fromrect(x, y, z):
    ra = numexpr.evaluate('arctan2(y,x)*57.295779513082323')
    dec = numexpr.evaluate('57.295779513082323*arctan2(z,sqrt(x**2+y**2))')
    return ra, dec


def rotation_matrix(rapol, decpol, ra0):
    """
    Return the rotation matrix corresponding to the pole of rapol, decpol
    and with zero of new latitude corresponding to ra=ra0
    This matrix need to be np.dot'ed with the input vector to get
    forward transform
    """
    tmppol = np.array(torect(rapol, decpol))  # pole axis
    tmpvec1 = np.array(torect(ra0, 0))  # x axis
    tmpvec1 = np.array(tmpvec1)

    tmpvec1[2] = (-tmppol[0] * tmpvec1[0] - tmppol[1] * tmpvec1[1]) / tmppol[2]
    tmpvec1 /= np.sqrt((tmpvec1**2).sum())
    tmpvec2 = np.cross(tmppol, tmpvec1)  # y axis
    M = np.array([tmpvec1, tmpvec2, tmppol])
    return M


def sphere_rotate(ra,
                  dec,
                  rapol=None,
                  decpol=None,
                  ra0=None,
                  revert=False,
                  mat=None):
    """ rotate ra,dec to a new spherical coordinate system where the pole is
    at rapol,decpol and the zeropoint is at ra=ra0
    revert flag allows to reverse the transformation
    """

    x, y, z = torect(ra, dec)
    if rapol is not None and decpol is not None and ra0 is not None:
        M = rotation_matrix(rapol, decpol, ra0)
    else:
        if mat is None:
            raise Exception('matrix or rapol,decpol, ra0'
                            ' need to be provided')
        else:
            M = mat

    if not revert:
        Axx, Axy, Axz = M[0]
        Ayx, Ayy, Ayz = M[1]
        Azx, Azy, Azz = M[2]
    else:
        Axx, Ayx, Azx = M[0]
        Axy, Ayy, Azy = M[1]
        Axz, Ayz, Azz = M[2]
    xnew = numexpr.evaluate('x*Axx+y*Axy+z*Axz')
    ynew = numexpr.evaluate('x*Ayx+y*Ayy+z*Ayz')
    znew = numexpr.evaluate('x*Azx+y*Azy+z*Azz')
    del x, y, z
    tmp = fromrect(xnew, ynew, znew)
    return (tmp[0], tmp[1])

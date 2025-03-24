import numpy as np
import matplotlib.pyplot as plt
import healpy

flip = -1  # for astronomical convention


class Formatter:

    def __init__(self, ra_shift=0):
        self.ra_shift = ra_shift

    def __call__(self, lon, lat):
        lon, lat = np.rad2deg([lon, lat])
        return '%f,%f' % (flip * lon + self.ra_shift, lat)


class Info:

    def __init__(self, ra_shift):
        self.ra_shift = ra_shift


def make_axes(ra_shift=None, dra_lab=30):
    """ create the axes
    dra_lab: is the step in ra labels
    ra_shift: this is sphere rotation in ra
    """
    if ra_shift is None:
        ra_shift = 0
    ax = plt.subplot(111, projection='mollweide')
    xticks = np.arange(0, 360, dra_lab)
    xticks_shift = (flip * xticks + ra_shift + 360 + 180) % 360 - 180
    xtick_labels = xticks.astype(int)
    xtick_labels = ['{:d}$^o$'.format(i) for i in xtick_labels]
    plt.xticks(ticks=np.radians(xticks_shift), labels=xtick_labels, rotation=0)
    plt.grid(True)
    ax.format_coord = Formatter(ra_shift)
    print(ra_shift)
    ax._extra = Info(ra_shift)


def scatter(ra, dec, ra_shift=None, dra_lab=30, overplot=False, **kwargs):
    """
    Scatter plot in Mollweide
    Parameters
    ----------
    ra : array
        right ascension in deg
    dec : array
        declination in deg
    ra_shift: shift/rotation in right ascension
    """
    if not overplot:
        make_axes(ra_shift=ra_shift, dra_lab=dra_lab)
    else:
        ra_shift = plt.gca()._extra.ra_shift
    ra_rad = np.deg2rad(flip * np.asarray(ra) + ra_shift)
    ra_rad = (ra_rad + 11 * np.pi) % (2 * np.pi) - np.pi
    dec_rad = np.deg2rad(dec)

    plt.scatter(ra_rad, dec_rad, **kwargs)


def _line_unwrapper(ra, dec):
    """
    when the line crosses the lon=-pi or lon=pi lines
    put nan there to avoid wrapping of lines
    """
    dra = np.diff(ra)
    cross = np.abs(dra) > np.pi
    x1 = _insert_nans(ra, cross)
    y1 = _insert_nans(dec, cross)
    return x1, y1


def _insert_nans(x, ind):
    """
    x   : original 1D array of length N
    ind : boolean array of length N-1
    """
    N = len(x)
    # Number of True entries tells us how many NaNs to insert
    num_nan = np.count_nonzero(ind)

    # The new array will have length N + num_nan
    x1 = np.full(N + num_nan, np.nan, dtype=float)

    # Build a "shift array" that tells us how many NaNs
    # should have been inserted up to each index in x
    # e.g. if ind = [True, False], cumsum(ind) = [1, 1]
    # we prepend 0 so we can handle indexing nicely
    shift_array = np.insert(np.cumsum(ind), 0, 0)

    # The position where x[i] should go in x1 is i + shift_array[i]
    positions = np.arange(N) + shift_array

    # Place the original x values into the correct positions
    x1[positions] = x

    return x1


def plot(ra, dec, ra_shift=None, dra_lab=30, overplot=False, **kwargs):
    """
    Plot/Line plot in Mollweide
    Parameters
    ----------
    ra : array
        right ascension in deg
    dec : array
        declination in deg
    ra_shift: shift in right ascension
    """
    if not overplot:
        make_axes(ra_shift=ra_shift, dra_lab=dra_lab)
        ra_shift = ra_shift or 0
    else:
        ra_shift = plt.gca()._extra.ra_shift
    ra_rad = np.deg2rad(flip * np.asarray(ra) + ra_shift)
    ra_rad = (ra_rad + 11 * np.pi) % (2 * np.pi) - np.pi
    dec_rad = np.deg2rad(dec)
    ra_rad1, dec_rad1 = _line_unwrapper(ra_rad, dec_rad)
    plt.plot(ra_rad1, dec_rad1, **kwargs)


def hpx_show(im,
             ra_shift=None,
             overplot=False,
             nest=True,
             pix_per_deg=10,
             **kwargs):
    """
    Plot a HEALPIX array on the mollweide map
    im: ndarray
        Array that needs to be plotted (with the length 12*nside**2)
    ra_shift: rotation in right ascension
    overplot: bool
        To reuse the axes or not
    pix_per_deg: float
        How fine should be the discretisation of the healpix map
    """
    if not overplot:
        make_axes(ra_shift=ra_shift)
        ra_shift = ra_shift or 0
    else:
        ra_shift = plt.gca()._extra.ra_shift
    nlon = int(360 * pix_per_deg)
    nlat = int(180 * pix_per_deg)
    lon_rad = np.radians(np.linspace(-180, 180, nlon))
    lat_rad = np.radians(np.linspace(-90, 90, nlat))

    lon2d, lat2d = np.meshgrid(lon_rad, lat_rad)
    nside = int(round((len(im) / 12)**.5))
    assert len(im) == 12 * nside**2
    hpx = healpy.ang2pix(nside,
                         flip * np.rad2deg(lon2d) + ra_shift,
                         np.rad2deg(lat2d),
                         lonlat=True,
                         nest=True)
    arr = im[hpx]
    ax = plt.gca()
    R = ax.pcolormesh(lon2d, lat2d, arr, **kwargs)
    return R

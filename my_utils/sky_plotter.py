import numpy as np
import matplotlib.pyplot as plt
import healpy

# Flip factor for astronomical convention (flips the longitude direction)
flip = -1


class Formatter:
    """
    Formatter for displaying coordinates in the Mollweide projection.
    Converts longitude and latitude from radians to degrees and applies
    a shift in right ascension (RA).
    """

    def __init__(self, ra_shift=0):
        """
        Parameters
        ----------
        ra_shift : float
            Shift in right ascension (degrees).
        """
        self.ra_shift = ra_shift

    def __call__(self, lon, lat):
        """
        Format the coordinates for display.

        Parameters
        ----------
        lon : float
            Longitude in radians.
        lat : float
            Latitude in radians.

        Returns
        -------
        str
            Formatted string with longitude and latitude in degrees.
        """
        lon, lat = np.rad2deg([lon, lat])
        return '%f,%f' % (flip * lon + self.ra_shift, lat)


class Info:
    """
    Container for storing additional information about the plot,
    such as the RA shift.
    """

    def __init__(self, ra_shift):
        """
        Parameters
        ----------
        ra_shift : float
            Shift in right ascension (degrees).
        """
        self.ra_shift = ra_shift


def make_axes(ra_shift=None, dra_lab=30):
    """
    Create a Mollweide projection plot with customized RA labels.

    Parameters
    ----------
    ra_shift : float, optional
        Shift in right ascension (degrees). Default is 0.
    dra_lab : int, optional
        Step size for RA labels (degrees). Default is 30.
    """
    if ra_shift is None:
        ra_shift = 0
    ax = plt.subplot(111, projection='mollweide')

    # Generate RA tick positions and labels
    xticks = np.arange(0, 360, dra_lab)
    xticks_shift = (flip * xticks + ra_shift + 360 + 180) % 360 - 180
    xtick_labels = xticks.astype(int)
    xtick_labels = ['{:d}$^o$'.format(i) for i in xtick_labels]

    # Set RA ticks and labels
    plt.xticks(ticks=np.radians(xticks_shift), labels=xtick_labels, rotation=0)
    plt.grid(True)

    # Set custom coordinate formatter
    ax.format_coord = Formatter(ra_shift)

    # Store additional information in the axes object
    ax._extra = Info(ra_shift)


def _line_unwrapper(ra, dec):
    """
    Handle line wrapping at the edges of the Mollweide projection.

    Parameters
    ----------
    ra : array-like
        Right ascension (radians).
    dec : array-like
        Declination (radians).

    Returns
    -------
    tuple
        Arrays of RA and Dec with NaNs inserted to avoid wrapping.
    """
    dra = np.diff(ra)
    cross = np.abs(dra) > np.pi  # Detect line crossings
    x1 = _insert_nans(ra, cross)
    y1 = _insert_nans(dec, cross)
    return x1, y1


def _insert_nans(x, ind):
    """
    Insert NaNs into an array to handle discontinuities.

    Parameters
    ----------
    x : array-like
        Original 1D array.
    ind : array-like
        Boolean array indicating where to insert NaNs.

    Returns
    -------
    array
        New array with NaNs inserted.
    """
    N = len(x)
    num_nan = np.count_nonzero(ind)  # Number of NaNs to insert

    # Create a new array with space for NaNs
    x1 = np.full(N + num_nan, np.nan, dtype=float)

    # Compute positions for original values
    shift_array = np.insert(np.cumsum(ind), 0, 0)
    positions = np.arange(N) + shift_array

    # Insert original values into the new array
    x1[positions] = x

    return x1


def scatter(ra, dec, ra_shift=None, dra_lab=30, overplot=False, **kwargs):
    """
    Scatter plot in Mollweide projection.

    Parameters
    ----------
    ra : array-like
        Right ascension (degrees).
    dec : array-like
        Declination (degrees).
    ra_shift : float, optional
        Shift in right ascension (degrees). Default is None.
    dra_lab : int, optional
        Step size for RA labels (degrees). Default is 30.
    overplot : bool, optional
        If True, reuse the existing axes. Default is False.
    kwargs : dict
        Additional arguments passed to `plt.scatter`.
    """
    if not overplot:
        make_axes(ra_shift=ra_shift, dra_lab=dra_lab)
    else:
        ra_shift = plt.gca()._extra.ra_shift

    # Convert RA and Dec to radians and apply RA shift
    ra_rad = np.deg2rad(flip * np.asarray(ra) + ra_shift)
    ra_rad = (ra_rad + 11 * np.pi) % (2 *
                                      np.pi) - np.pi  # Wrap RA to [-pi, pi]
    dec_rad = np.deg2rad(dec)

    # Plot the scatter points
    plt.scatter(ra_rad, dec_rad, **kwargs)


def plot(ra, dec, ra_shift=None, dra_lab=30, overplot=False, **kwargs):
    """
    Line plot in Mollweide projection.

    Parameters
    ----------
    ra : array-like
        Right ascension (degrees).
    dec : array-like
        Declination (degrees).
    ra_shift : float, optional
        Shift in right ascension (degrees). Default is None.
    dra_lab : int, optional
        Step size for RA labels (degrees). Default is 30.
    overplot : bool, optional
        If True, reuse the existing axes. Default is False.
    kwargs : dict
        Additional arguments passed to `plt.plot`.
    """
    if not overplot:
        make_axes(ra_shift=ra_shift, dra_lab=dra_lab)
        ra_shift = ra_shift or 0
    else:
        ra_shift = plt.gca()._extra.ra_shift

    # Convert RA and Dec to radians and apply RA shift
    ra_rad = np.deg2rad(flip * np.asarray(ra) + ra_shift)
    ra_rad = (ra_rad + 11 * np.pi) % (2 *
                                      np.pi) - np.pi  # Wrap RA to [-pi, pi]
    dec_rad = np.deg2rad(dec)

    # Handle line wrapping
    ra_rad1, dec_rad1 = _line_unwrapper(ra_rad, dec_rad)

    # Plot the line
    plt.plot(ra_rad1, dec_rad1, **kwargs)


def hpx_show(im,
             ra_shift=None,
             overplot=False,
             nest=True,
             pix_per_deg=10,
             **kwargs):
    """
    Plot a HEALPix array on a Mollweide projection.

    Parameters
    ----------
    im : ndarray
        HEALPix array (length 12 * nside^2).
    ra_shift : float, optional
        Shift in right ascension (degrees). Default is None.
    overplot : bool, optional
        If True, reuse the existing axes. Default is False.
    nest : bool, optional
        If True, use nested HEALPix ordering. Default is True.
    pix_per_deg : float, optional
        Resolution of the plot (pixels per degree). Default is 10.
    kwargs : dict
        Additional arguments passed to `pcolormesh`.

    Returns
    -------
    QuadMesh
        The plotted HEALPix map.
    """
    if not overplot:
        make_axes(ra_shift=ra_shift)
        ra_shift = ra_shift or 0
    else:
        ra_shift = plt.gca()._extra.ra_shift

    # Generate a grid of longitude and latitude
    nlon = int(360 * pix_per_deg)
    nlat = int(180 * pix_per_deg)
    lon_rad = np.radians(np.linspace(-180, 180, nlon))
    lat_rad = np.radians(np.linspace(-90, 90, nlat))
    lon2d, lat2d = np.meshgrid(lon_rad, lat_rad)

    # Compute HEALPix pixel indices
    nside = int(round((len(im) / 12)**0.5))
    assert len(im) == 12 * nside**2
    hpx = healpy.ang2pix(nside,
                         flip * np.rad2deg(lon2d) + ra_shift,
                         np.rad2deg(lat2d),
                         lonlat=True,
                         nest=nest)

    # Map HEALPix values to the grid
    arr = im[hpx]

    # Plot the HEALPix map
    ax = plt.gca()
    R = ax.pcolormesh(lon2d, lat2d, arr, **kwargs)
    return R

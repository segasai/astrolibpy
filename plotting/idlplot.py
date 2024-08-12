# Copyright (C) 2009-2024 Sergey Koposov
# This file is part of astrolibpy
#
# astrolibpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# astrolibpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with astrolibpy.  If not, see <http://www.gnu.org/licenses/>.

import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.stats
import scipy.ndimage.filters
import scipy.stats
import matplotlib
import math
import copy
import warnings

# this module is by default in interactive regime
plt.ion()


def listToArr(x):
    if isinstance(x, np.ma.MaskedArray):
        x = x.filled(np.nan)
    return np.asarray(x)


def listToArrFlat(x):
    if isinstance(x, np.ma.MaskedArray):
        x = x.filled(np.nan)
    return np.asarray(x).flatten()


def _get_marker(ps, linestyle):
    """
    Wrapper for point markers which understand idl-like ps options
    (e.g. ps=3 for the point, ps=4 for the diamond)
    """
    outlinestyle = ' '
    markerHash = {3: '.', 4: 'D', 7: '+', 2: '*', 6: 's'}
    try:
        marker = markerHash[ps]
    except KeyError:
        if isinstance(ps, str):
            marker = ps
        else:
            marker = ' '
            if linestyle is not None:
                outlinestyle = linestyle
            else:
                outlinestyle = '-'
    return (marker, outlinestyle)


def _get_ranges(x, y, xr, yr, xrange, yrange, pad_range):
    if xr is not None and xrange is not None:
        raise ValueError('xr and xrange are both specified')
    if yr is not None and yrange is not None:
        raise ValueError('yr and yrange are both specified')
    xrange = xrange or xr
    yrange = yrange or yr

    if xrange is None and yrange is None:
        ind = np.isfinite(x)
        if not ind.any():
            xrange = [0, 1]
        else:
            xrange = [np.min(x[ind]), np.max(x[ind])]
        ind = np.isfinite(y)
        if not ind.any():
            yrange = [0, 1]
        else:
            yrange = [np.min(y[ind]), np.max(y[ind])]
        assert (pad_range >= 0)
        xrange = [
            xrange[0] - pad_range * (xrange[1] - xrange[0]),
            xrange[1] + pad_range * (xrange[1] - xrange[0])
        ]
        yrange = [
            yrange[0] - pad_range * (yrange[1] - yrange[0]),
            yrange[1] + pad_range * (yrange[1] - yrange[0])
        ]

        del ind
    elif xrange is None and yrange is not None:
        ind = (y < max(yrange[1], yrange[0])) & (y > min(
            yrange[0], yrange[1])) & np.isfinite(x)
        if ind.any():
            xrange = [np.min(x[ind]), np.max(x[ind])]
        else:
            xrange = [np.min(x), np.max(x)]
        del ind
    elif xrange is not None and yrange is None:
        ind = (x < np.maximum(xrange[1], xrange[0])) & (x > np.minimum(
            xrange[0], xrange[1])) & np.isfinite(y)
        if ind.any():
            yrange = [np.min(y[ind]), np.max(y[ind])]
        else:
            yrange = [np.min(y), np.max(y)]
        del ind
    if len(yrange) != 2 or len(xrange) != 2:
        raise ValueError("Wrong xrange or yrange")
    return xrange, yrange


def filter_epa(im, kernsize):
    if not hasattr(kernsize, '__iter__'):
        kernsize1, kernsize2 = kernsize, kernsize
    else:
        kernsize1, kernsize2 = kernsize[0:2]
    k11 = math.ceil(kernsize1)
    k12 = math.ceil(kernsize2)
    xgrid, ygrid = np.mgrid[-k11:k11:1, -k12:k12:1]
    r2 = (xgrid / kernsize1)**2 + (ygrid / kernsize2)**2
    filt = (1 - r2) * (r2 <= 1)
    filt = filt / filt.sum()
    # IMPORTANT TRANSPOSITION, because the image first dimension is y
    im1 = scipy.ndimage.filters.convolve(im, filt.T, mode='reflect')
    return im1


def smoother(arr, smooth=None, kernel=None):
    if smooth is not None:
        if kernel == 'gau':
            if hasattr(smooth, '__iter__'):
                smooth = smooth[::-1]
                # because gaussian_filter convention is second dimension is x
            arr = scipy.ndimage.filters.gaussian_filter(arr * 1., smooth)
        elif kernel == 'epa':
            arr = filter_epa(arr * 1., smooth)
        else:
            raise Exception('Wrong kernel')
    return arr


def exceptionDecorator(func):
    """
    The purpose of this decorator is to
    do the plotting in the non-interactive mode
    and then return back to whatever previous mode was
    (irrespective if there was an exception or not)
    """

    def wrapper(*args, **kwargs):
        with plt.ioff():
            ret = func(*args, **kwargs)
        if plt.isinteractive():
            plt.gcf().show()
        return ret

    wrapper.__doc__ = func.__doc__

    return wrapper


def __findKnuth(x, minv, maxv):
    """
    Implement Knuth method for histogram bin selection

    """
    N = ((x >= minv) & (x <= maxv)).sum()

    def funcer(M):
        hh, loc = np.histogram(x, bins=M, range=[minv, maxv])
        return np.log(M) + 1. / N * (scipy.special.gammaln(M / 2.) -
                                     M * scipy.special.gammaln(0.5) -
                                     scipy.special.gammaln(N + M / 2.) +
                                     scipy.special.gammaln(hh + 0.5).sum())

    maxN = 1000
    ns = np.arange(1, maxN + 1)
    vals = ns * 0.
    for i in range(len(ns)):
        vals[i] = funcer(ns[i])
    bestn = ns[np.argmax(vals)]
    if bestn == maxN:
        print('WARNING the best number of bins is > maxbin(%d)' % (maxN))
    return bestn


@exceptionDecorator
def plothist(x,
             bin=None,
             nbins=None,
             xrange=None,
             yrange=None,
             min=None,
             max=None,
             overplot=False,
             color='black',
             xlog=False,
             ylog=False,
             nan=False,
             weights=None,
             norm=False,
             kernel=None,
             retpoints=False,
             adaptive=False,
             adaptive_thresh=30,
             adaptive_depth=[2, 10],
             weight_norm=False,
             apply_func=None,
             statistic=None,
             knuth=False,
             cumulative=False,
             **kw):
    """
    Plot the 1D histogram
    Example:
    >> plothist(dat, bin=0.1, min=0, max=3)

    Keyword parameters:
    ------------------
    bin
        the binsize(float)
    nbins
        number of bins(integer)
        It cannot be specified together with the bin= parameter
    xlog, ylog
        log the appropriate axis
    weights
        the 1-D array of weights used in the histogram creation
    nan
        boolean flag to ignore nan's
    norm
        boolean flag to normalize the histogram by the peak value
    min,max
        range of data for which the histogram is constructed
    retpoints
        boolean parameter controlling whether to return or not the
        computed histogram. If yes the tuple with two arrays
        (bin centers, Number of points in bins) is returned
    overplot
        boolean parameter for overplotting
    adaptive
        boolean for turning on/off the adaptive regime of
        histogramming (adaptive bin size).
        If True weights, nbins, bin,kernel parameters are ignored
    adaptive_thresh
        the limiting number of points in the bin for the adaptive
        histogramming (default 30)
    adaptive_depth
        the list of two integers for the detalisation levels of
        adaptive histogramming (default [2,10])
    weight_norm
        if True the value in each bin is mean weight of points within
        the bin
    """
    if nan:
        ind = np.isfinite(x)
        if weights is not None:
            ind = np.isfinite(weights) & ind
        dat = x[ind]
        if weights is not None:
            weights = weights[ind]
    else:
        dat = x
    maxNone = False
    if min is None:
        min = np.nanmin(dat)
    if max is None:
        maxNone = True
        max = np.nanmax(dat)

    if bin is None and nbins is None:
        if not knuth:
            nbins = 100
        else:
            nbins = __findKnuth(dat, min, max)
        bin = (max - min) * 1. / nbins
    elif nbins is None:
        nbins = int(math.ceil((max - min) * 1. / bin))
        if maxNone:
            max = min + nbins * bin
    elif bin is None:
        bin = (max - min) * 1. / nbins
    else:
        warnings.warn(
            'both bin= and nbins= keywords were specified in ' +
            'the plothist call', RuntimeWarning)
        pass
        # if both nbins and bin are defined I don't do anything
        # it may be non-intuitive if kernel option is used, because
        # it uses both nbins and bin options
    if cumulative:
        if (kernel is not None or adaptive or weights is not None):
            raise RuntimeError(
                'cumulative is incompatible with weights, kernel ' +
                'or adaptive options')
    if kernel is None:
        if not adaptive:
            if np.isscalar(weights) and weights is not None:
                weights = np.zeros_like(dat) + weights
            if statistic is None:
                hh, loc = np.histogram(dat,
                                       range=(min, max),
                                       bins=nbins,
                                       weights=weights)

                if weight_norm:
                    hh1, loc = np.histogram(dat,
                                            range=(min, max),
                                            bins=nbins,
                                            weights=None)
                    hh = hh * 1. / hh1
            else:
                S = scipy.stats.binned_statistic(dat,
                                                 weights,
                                                 statistic,
                                                 range=[min, max],
                                                 bins=nbins)
                hh = S.statistic
                loc = S.bin_edges
        else:
            import adabinner
            hh, loc = adabinner.hist(dat,
                                     xmin=min,
                                     xmax=max,
                                     hi=adaptive_depth,
                                     thresh=adaptive_thresh)
        if cumulative:
            hh = np.cumsum(hh)
        hh1 = np.repeat(hh, 2)
        loc1 = np.concatenate(([loc[0]], np.repeat(loc[1:-1], 2), [loc[-1]]))
    else:
        loc1 = np.linspace(min, max, nbins * 10)
        import sklearn.neighbors
        xind = np.isfinite(dat)
        kde = sklearn.neighbors.KernelDensity(bandwidth=bin, kernel=kernel)
        kde.fit(np.asarray(dat[xind]).flatten().reshape(-1, 1))
        hh1 = np.exp(kde.score_samples(loc1.reshape(-1, 1)))
        if weights is not None:
            print('WARNING weights ignored for KDE !')
    if overplot:
        func = oplot
    else:
        func = plot
    if norm:
        hh1 = hh1 * 1. / hh1.max()
    kw['ps'] = kw.get('ps') or 0
    if 'yr' not in kw:
        kw['yr'] = [np.nanmin(hh1), np.nanmax(hh1)]
    if 'xr' not in kw:
        kw['xr'] = [min, max]
    if apply_func is not None:
        hh1 = apply_func(loc1, hh1)
    func(loc1, hh1, color=color, xlog=xlog, ylog=ylog, **kw)
    if retpoints:
        return 0.5 * (loc[1:] + loc[:-1]), hh


@exceptionDecorator
def plot(arg1,
         arg2=None,
         xrange=None,
         yrange=None,
         ps=0,
         thick=None,
         xtitle=None,
         ytitle=None,
         color='black',
         noerase=False,
         overplot=False,
         position=None,
         ylog=False,
         xlog=False,
         xr=None,
         yr=None,
         title=None,
         label=None,
         nodata=False,
         linestyle=None,
         markersize=None,
         xaxis_formatter=None,
         yaxis_formatter=None,
         autoscalex=False,
         autoscaley=False,
         markerfacecolor=None,
         markeredgecolor=None,
         markeredgewidth=None,
         axis=None,
         pad_range=0,
         **kw):
    """ Plot your data in an IDL-like way with a lot of options in one
    command.
    Example:

    >> plot(x,y,xrange=[0,39],yrange=[-1,10],ps=4,xtitle="X",
        color='black',position=[0.1,0.1,0.9,0.9], xlog=True)
    """

    if arg2 is None:
        y = listToArrFlat(arg1)
        x = np.arange(len(y))
    else:
        x = listToArrFlat(arg1)
        y = listToArrFlat(arg2)

    if not noerase:
        plt.gcf().clf()
    if position is not None and axis is None:
        mypos = position[:]
        mypos[2] = position[2] - position[0]
        mypos[3] = position[3] - position[1]
        plt.axes(mypos)
    if axis is None:
        axis = plt.gca()

    marker, linestyle = _get_marker(ps, linestyle)
    xrange, yrange = _get_ranges(x, y, xr, yr, xrange, yrange, pad_range)

    if not overplot:
        axis.minorticks_on()
    if xtitle is not None:
        axis.set_xlabel(xtitle)
    if ytitle is not None:
        axis.set_ylabel(ytitle)

    axis.set_autoscalex_on(autoscalex)
    axis.set_autoscaley_on(autoscaley)
    if not overplot:
        axis.set_xlim(xrange)
        axis.set_ylim(yrange)

    if xlog:
        axis.set_xscale('log', subs=[2, 3, 4, 5, 6, 7, 8, 9])
    if ylog:
        axis.set_yscale('log', subs=[2, 3, 4, 5, 6, 7, 8, 9])
    if xaxis_formatter is not None:
        axis.xaxis.set_major_formatter(xaxis_formatter)
    if yaxis_formatter is not None:
        axis.yaxis.set_major_formatter(yaxis_formatter)

    if title is not None:
        plt.title(title)
    if not nodata:
        axis.plot(x,
                  y,
                  marker=marker,
                  linestyle=linestyle,
                  linewidth=thick,
                  color=color,
                  label=label,
                  markersize=markersize,
                  markerfacecolor=markerfacecolor,
                  markeredgecolor=markeredgecolor,
                  markeredgewidth=markeredgewidth,
                  **kw)


@exceptionDecorator
def plot_scatter(x,
                 y,
                 c=None,
                 s=None,
                 xrange=None,
                 yrange=None,
                 thick=None,
                 xtitle=None,
                 ytitle=None,
                 color=None,
                 noerase=False,
                 overplot=False,
                 position=None,
                 ylog=False,
                 xlog=False,
                 xr=None,
                 yr=None,
                 title=None,
                 label=None,
                 marker='o',
                 markersize=None,
                 xaxis_formatter=None,
                 yaxis_formatter=None,
                 autoscalex=False,
                 autoscaley=False,
                 axis=None,
                 bar=False,
                 bar_label=None,
                 bar_fraction=0.05,
                 pad_range=0,
                 **kw):
    """
    Plot your data in an IDL-like way with a lot of options in one
    command. This particular command is a wrapper around matplotlib.scatter
    Thus plot scatter plot where the colour and size of the points is
    determined by an array.
    Example:

    >> plot_scatter(x, y, c=z, xrange=[0,39],yrange=[-1,10],ps=4,xtitle="X",
        color='black',position=[0.1,0.1,0.9,0.9], xlog=True)

    Arguments:
    x, y: arrays to be plotted
    c: array of colours (can be None)
    s: array of sizes (can be None)
    """

    x = listToArrFlat(x)
    y = listToArrFlat(y)

    if c is not None:
        c = listToArrFlat(c)
    if s is not None:
        s = listToArrFlat(s)

    if not noerase:
        plt.gcf().clf()
    if position is not None and axis is None:
        mypos = position[:]
        mypos[2] = position[2] - position[0]
        mypos[3] = position[3] - position[1]
        plt.axes(mypos)
    if axis is None:
        axis = plt.gca()

    xrange, yrange = _get_ranges(x, y, xr, yr, xrange, yrange, pad_range)
    if not overplot:
        axis.minorticks_on()
    if xtitle is not None:
        axis.set_xlabel(xtitle)
    if ytitle is not None:
        axis.set_ylabel(ytitle)

    axis.set_autoscalex_on(autoscalex)
    axis.set_autoscaley_on(autoscaley)
    if not overplot:
        axis.set_xlim(xrange)
        axis.set_ylim(yrange)

    if xlog:
        axis.set_xscale('log', subs=[2, 3, 4, 5, 6, 7, 8, 9])
    if ylog:
        axis.set_yscale('log', subs=[2, 3, 4, 5, 6, 7, 8, 9])

    if xaxis_formatter is not None:
        axis.xaxis.set_major_formatter(xaxis_formatter)
    if yaxis_formatter is not None:
        axis.yaxis.set_major_formatter(yaxis_formatter)

    if title is not None:
        plt.title(title)
    scat = axis.scatter(x,
                        y,
                        c=c,
                        s=s,
                        marker=marker,
                        label=label,
                        color=color,
                        **kw)
    if bar:
        if int(''.join((matplotlib.__version__).split('.')[:2])) >= 11:
            kw = {'use_gridspec': True}
        else:
            kw = {}
        cb = plt.colorbar(fraction=bar_fraction, mappable=scat, **kw)
        if bar_label is not None:
            cb.set_label(bar_label)


def oplot(x, y=None, **kw):
    """
    Overplot your data (just a simple wrapper around plot)

    Example:
    >> oplot(x,2+y/10.,ps=3,color='blue')
    """

    plot(x, y, noerase=True, overplot=True, **kw)


@exceptionDecorator
def ploterror(x,
              y,
              err0,
              err1=None,
              color='black',
              ps=0,
              ecolor='black',
              overplot=False,
              noerase=False,
              elinewidth=None,
              capsize=None,
              markersize=None,
              markeredgecolor=None,
              markerfacecolor=None,
              autoscalex=False,
              autoscaley=False,
              label=None,
              **kw):
    """
    Plot the data with error-bars

    Example:
    >> ploterror(x,y,erry) # plot only Y error-bars
    >> ploterror(x,y,errx,erry) # plot data with X and Y error-bars

    Keyword parameters:
    ------------------
    capsize
        integer param controlling the size of hats/caps of the error-bars
    ecolor
        color of the error-bars (different from the main color)
    """
    if overplot:
        noerase = True
    if err1 is None:
        erry = listToArr(err0)
    else:
        erry = listToArr(err1)
        errx = listToArr(err0)
    x = listToArrFlat(x)
    y = listToArrFlat(y)
    kw0 = kw.copy()
    if kw0.get('yr') is None:
        kw0['yr'] = [np.nanmin(y - erry), np.nanmax(y + erry)]

    if markersize is not None:
        kw0['markersize'] = markersize
    if markeredgecolor is not None:
        kw0['markeredgecolor'] = markeredgecolor
    if markerfacecolor is not None:
        kw0['markerfacecolor'] = markerfacecolor
    plot(x, y, color=color, ps=ps, overplot=overplot, noerase=noerase, **kw0)
    (marker, outlinestyle) = _get_marker(ps, None)
    kw1 = {
        'ecolor': ecolor,
        'marker': marker,
        'color': color,
        'linestyle': outlinestyle,
        'elinewidth': elinewidth,
        'markeredgewidth': elinewidth
    }
    if 'alpha' in kw0:
        kw1['alpha'] = kw0['alpha']
    if 'zorder' in kw0:
        kw1['zorder'] = kw0['zorder']
    if label is not None:
        kw1['label'] = label

    if markersize is not None:
        kw1['markersize'] = markersize
    if markeredgecolor is not None:
        kw1['markeredgecolor'] = markeredgecolor
    if markerfacecolor is not None:
        kw1['markerfacecolor'] = markerfacecolor

    if capsize is not None:
        kw1['capsize'] = capsize

    plt.gca().set_autoscalex_on(autoscalex)
    plt.gca().set_autoscaley_on(autoscaley)

    if err1 is None:
        plt.gca().errorbar(x, y, erry, **kw1)
    else:
        plt.gca().errorbar(x, y, xerr=errx, yerr=erry, **kw1)


@exceptionDecorator
def tvaxis(image,
           xmin=None,
           xmax=None,
           ymin=None,
           ymax=None,
           xtitle="",
           ytitle="",
           title="",
           vmin=None,
           vmax=None,
           aspect="auto",
           xlog=False,
           ylog=False,
           position=None,
           noerase=False,
           bar=False,
           bar_label=None,
           bar_fraction=0.05,
           zlog=False,
           zsqrt=False,
           smooth=None,
           vmaxfrac=None,
           xflip=False,
           yflip=False,
           kernel='gau',
           vminfrac=None,
           **kw):
    """
    Display the 2D image with proper axes (similar to plt.imshow)
    Example:
    >> tvaxis(im,-20,10,-40,50)

    Keyword parameters:
    ------------------
    xmin,xmax,ymin,ymax
        the ranges for x,y where the histogram is constructed.
        These params can be specified not only as keywords but also as normal
        arguments
    vmin,vmax
        the ranges of intensities shown on the plot
    smooth
        if not None this parameter controls additional smoothing of the 2D
        histogram
    bar
        boolean parameter for switching on/off the plotting of the color-bar

    """

    if not noerase:
        plt.clf()

    if position is not None:
        mypos = position[:]
        mypos[2] = position[2] - position[0]
        mypos[3] = position[3] - position[1]
        plt.axes(mypos)
    if xmin is None:
        xmin = 0
    if ymin is None:
        ymin = 0
    if xmax is None:
        xmax = image.shape[0]
    if ymax is None:
        ymax = image.shape[1]
    if image.ndim == 3:
        im = np.transpose(image, axes=(1, 0, 2))
    elif image.ndim == 2:
        im = image.T
    else:
        raise ValueError('Wrong dimensions of the input array')
    im = smoother(im, smooth=smooth, kernel=kernel)

    vmin, vmax, norm = _set_vminmax_norm(vmin, vminfrac, vmax, vmaxfrac, zlog,
                                         zsqrt, im)

    if xflip:
        im = np.fliplr(im)
    if yflip:
        im = np.flipud(im)

    axim = plt.imshow(im,
                      extent=(xmin, xmax, ymin, ymax),
                      aspect=aspect,
                      norm=norm,
                      origin='lower',
                      **kw)
    if xlog:
        plt.gca().set_xscale('log')
    if ylog:
        plt.gca().set_yscale('log')

    plt.gca().set_xlabel(xtitle)
    plt.gca().set_ylabel(ytitle)
    plt.gca().minorticks_on()

    if title is not None:
        plt.title(title)
    if bar:
        if int(''.join((matplotlib.__version__).split('.')[:2])) >= 11:
            kw = {'use_gridspec': True}
        else:
            kw = {}
        cb = plt.colorbar(fraction=bar_fraction, **kw)
        if bar_label is not None:
            cb.set_label(bar_label)

    return axim


def _set_vminmax_norm(vmin, vminfrac, vmax, vmaxfrac, zlog, zsqrt, im):
    if vminfrac is not None and vmin is None:
        vmin = scipy.stats.scoreatpercentile(im[np.isfinite(im)],
                                             100 * vminfrac)
    if vmaxfrac is not None and vmax is None:
        vmax = scipy.stats.scoreatpercentile(im[np.isfinite(im)],
                                             100 * vmaxfrac)
    if vmin is not None and vmax is not None and vmin >= vmax:
        warnings.warn("vmin is >= vmax... Resetting their values",
                      RuntimeWarning)
        vmin = None
        vmax = None

    if zlog:
        # here comes a special case
        # if vmin is None and actual minimum is X then it is undistinguishable from
        # 0 in log scale, so we set it to 0.9*X
        if vmin is None:
            try:
                vmin = np.min(im[im > 0])
                vmin = 0.9 * vmin
            except ValueError:
                pass
        norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
    elif zsqrt:
        norm = matplotlib.colors.PowerNorm(0.5, vmin=vmin, vmax=vmax)
    else:
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    return vmin, vmax, norm


def __parse_limit(val, x, oper=None):
    if isinstance(val, str):
        if val[0] == 'p':
            perc = float(val[1:])
            return scipy.stats.scoreatpercentile(x[np.isfinite(x)], perc)
        else:
            raise ValueError('unsupported limit format %s' % (val))
    if val is not None:
        return val
    else:
        return oper(x)


@exceptionDecorator
def tvhist2d(x,
             y,
             xmin=None,
             xmax=None,
             ymin=None,
             ymax=None,
             vmin=None,
             vmax=None,
             bins=[100, 100],
             xtitle="",
             ytitle="",
             title=None,
             noerase=False,
             weights=None,
             zlog=False,
             xflip=False,
             yflip=False,
             smooth=None,
             quick=True,
             cmap='gray_r',
             normx=None,
             normy=None,
             xlog=False,
             ylog=False,
             weight_norm=False,
             vminfrac=None,
             vmaxfrac=None,
             position=None,
             xaxis_formatter=None,
             yaxis_formatter=None,
             bar=False,
             bar_label='',
             bar_fraction=0.05,
             bar_pad=0.05,
             bar_ticks_locator=None,
             bar_formatter=None,
             apply_func=None,
             zsqrt=False,
             ret_hist=False,
             interpolation='nearest',
             scatter_thresh=None,
             scatter_opt={},
             subplot=None,
             kernel='gau',
             statistic=None,
             **kw):
    """ Plot the 2D histogram of the data
    Example:
    >> tvhist2d(xs,ys,bins=[30,30])

    >> tvhist2d(xs,ys,0,10,-1,2,bins=[30,30])

    Keyword arguments:
    -----------------
    xmin,xmax,ymin,ymax
        the ranges for x,y where the histogram is constructed.
        These params can be specified not only as keywords but also as
        normal arguments
        >> tvhist2d(xs,ys,xmin=0,ymin=10,ymin=-1,ymax=2,bins=[30,30])
        >> tvhist2d(xs,ys,0,10,-1,2,bins=[30,30])
        Also xmin,..ymax can be specified as percentiles as strings, i.e.
        >> tvhist2d(xs,ys,'p1','p99','p5','p95') for 1/99/5/95% respectively
    vmin,vmax
        the ranges of intensities shown on the plot
    vminfrac, vmaxfrac
        The percenitiles for vmin,vmax (i.e. 0,1 correspond to 0%,100%)
    bins
        the list of two integers specifying how many bins in x,y you want
    smooth
        if not None this parameter controls additional smoothing of the 2D
        histogram
    xflip, yflip
        boolean parameters allowing to flip x,y axes
    normx, normy
        params controlling the normalization of the histogram along X or Y axes
        if normx=='sum' then the normalization is done in such way that the sum
        is the same for each row/column
        if normx=='max' then the normalization is don in such way that
        the brightest pixel value will be the same for each row/column
    weight_norm
        if True the value in each bin is mean weight of points within
        the bin
    weights
        The array of weights for weighted histogram or for computing a
        statatistic
    bar
        boolean parameter for switching on/off the plotting of the color-bar
    bar_pad, bar_fraction
        padding and fraction for bar placement
    bar_ticks_locator
        locator for the tickmarks on the colorbar
    statistic
        The function or statistic applied to weights. It can be 'mean',
        'median' or arbitrary function.

    """

    x1 = listToArrFlat(x)
    y1 = listToArrFlat(y)

    xmin = __parse_limit(xmin, x1, np.nanmin)
    xmax = __parse_limit(xmax, x1, np.nanmax)
    ymin = __parse_limit(ymin, y1, np.nanmin)
    ymax = __parse_limit(ymax, y1, np.nanmax)

    plot_range = (xmin, xmax, ymin, ymax)
    hist_range = [[ymin, ymax], [xmin, xmax]]
    hist_range_ordered = []
    hist_range_reverse = []
    for i in range(2):
        left, right = hist_range[i]
        if left > right:
            hist_range_ordered.append([right, left])
            hist_range_reverse.append(True)
        else:
            hist_range_ordered.append([left, right])
            hist_range_reverse.append(False)

    binsRev = bins[::-1]
    if statistic is None:
        if not quick:
            hh, yedges, xedges = np.histogram2d(y1,
                                                x1,
                                                range=hist_range_ordered,
                                                bins=binsRev,
                                                weights=weights)
            if weight_norm:
                hh1, yedges, xedges = np.histogram2d(y1,
                                                     x1,
                                                     range=hist_range_ordered,
                                                     bins=binsRev,
                                                     weights=None)
                hh = hh * 1. / hh1

        else:
            import quick_hist
            hh = quick_hist.quick_hist((y1, x1),
                                       range=hist_range_ordered,
                                       nbins=binsRev,
                                       weights=weights)
            if weight_norm:
                hh1 = quick_hist.quick_hist((y1, x1),
                                            range=hist_range_ordered,
                                            nbins=binsRev)
                hh = hh * 1. / hh1
    else:
        hh = scipy.stats.binned_statistic_2d(y1,
                                             x1,
                                             weights,
                                             statistic,
                                             range=hist_range_ordered,
                                             bins=binsRev).statistic
    # reverse some of the dimensions if the axis range is reverse
    slices = []
    for i in range(2):
        if hist_range_reverse[i]:
            slices.append(slice(None, None, -1))
        else:
            slices.append(slice(None, None, None))
    hh = hh[tuple(slices)]
    if apply_func is not None:
        hh = apply_func(hh)
    hh = smoother(hh, smooth=smooth, kernel=kernel)

    if normx is not None:
        if normx == 'sum':
            hhs = hh.sum(axis=0)
            hhs = hhs + (hhs == 0)
            hh = hh * 1. / hhs[None, :]
        elif normx == 'max':
            hhs = np.max(hh, axis=0)
            hhs = hhs + (hhs == 0)
            hh = hh * 1. / hhs[None, :]
        elif normx == 'median':
            hhs = np.median(hh, axis=0)
            hhs = hhs + (hhs == 0)
            hh = hh * 1. / hhs[None, :]
        else:
            raise Exception('unknown normx mode')

    if normy is not None:
        if normy == 'sum':
            hhs = hh.sum(axis=1)
            hhs = hhs + (hhs == 0)
            hh = hh * 1. / hhs[:, None]
        elif normy == 'max':
            hhs = np.max(hh, axis=1)
            hhs = hhs + (hhs == 0)
            hh = hh * 1. / hhs[:, None]
        elif normy == 'median':
            hhs = np.median(hh, axis=1)
            hhs = hhs + (hhs == 0)
            hh = hh * 1. / hhs[:, None]
        else:
            raise Exception('unknown normy mode')

    if scatter_thresh is not None:
        locx = np.linspace(xmin, xmax, bins[0] + 1, True)
        locy = np.linspace(ymin, ymax, bins[1] + 1, True)
        posx = np.digitize(x1, locx)
        posy = np.digitize(y1, locy)
        # select points within the histogram
        ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
        hhsub = hh.T[posx[ind] - 1, posy[ind] -
                     1]  # values of the histogram where the points are
        x1 = x1[ind][hhsub < scatter_thresh]  # low density points
        y1 = y1[ind][hhsub < scatter_thresh]
        hh = hh.astype(np.float)
        hh[hh <
           scatter_thresh] = np.nan  # fill the areas with low density by NaNs

    if xflip:
        plot_range = tuple(plot_range[_] for _ in [1, 0, 2, 3])
        hh = np.fliplr(hh)
    if yflip:
        plot_range = tuple(plot_range[_] for _ in [0, 1, 3, 2])
        hh = np.flipud(hh)
    if subplot is not None:
        noerase = True
        if hasattr(subplot, '__iter__'):
            plt.subplot(*subplot)
        else:
            plt.subplot(subplot)
    if not noerase:
        plt.gcf().clf()
    if position is not None:
        mypos = position[:]
        mypos[2] = position[2] - position[0]
        mypos[3] = position[3] - position[1]
        plt.axes(mypos)
    axis = plt.gca()
    axis.set_xlabel(xtitle)
    axis.set_ylabel(ytitle)
    axis.minorticks_on()

    if xaxis_formatter is not None:
        axis.xaxis.set_major_formatter(xaxis_formatter)
    if yaxis_formatter is not None:
        axis.yaxis.set_major_formatter(yaxis_formatter)

    vmin, vmax, norm = _set_vminmax_norm(vmin, vminfrac, vmax, vmaxfrac, zlog,
                                         zsqrt, hh)

    axim = plt.imshow(hh,
                      extent=plot_range,
                      aspect='auto',
                      interpolation=interpolation,
                      cmap=cmap,
                      norm=norm,
                      origin='lower',
                      **kw)
    if scatter_thresh is not None:
        if 'ps' not in scatter_opt:
            scatter_opt = copy.copy(scatter_opt)
            scatter_opt['ps'] = 3
        oplot(x1, y1, **scatter_opt)

    if bar:
        if bar_formatter is None:
            bar_formatter = matplotlib.ticker.ScalarFormatter()
        cb = plt.colorbar(fraction=bar_fraction,
                          pad=bar_pad,
                          ax=axis,
                          format=bar_formatter,
                          ticks=bar_ticks_locator,
                          use_gridspec=True)
        if bar_label is not None:
            cb.set_label(bar_label)
    if xlog:
        plt.gca().set_xscale('log')
    if ylog:
        plt.gca().set_yscale('log')
    if title is not None:
        plt.title(title)
    if ret_hist:
        return hh

    return axim


@exceptionDecorator
def contour(z,
            x=None,
            y=None,
            xrange=None,
            yrange=None,
            zrange=None,
            xr=None,
            yr=None,
            zr=None,
            title=None,
            xtitle=None,
            ytitle=None,
            position=None,
            xlog=False,
            ylog=False,
            zlog=False,
            xticklabel=None,
            yticklabel=None,
            zticklabel=None,
            levels=None,
            nlevels=256,
            c_color="black",
            cm="jet",
            c_line="solid",
            c_levels=None,
            c_charsize=12.0,
            c_thick=1,
            thick=1,
            font="monospace",
            weight="normal",
            charsize=14.0,
            bar=True,
            fill=True,
            overplot=False,
            noerase=False,
            c_label=False,
            bar_fraction=0.05,
            xaxis_formatter=None,
            label=None,
            yaxis_formatter=None):
    """
    Plot the contours of the 2d array.
    Example:
    >> contour(z)
    if you have the x and y coordinates of your array then you can use them
    >> contour(z, x, y)
    In that case  x,y can be either 2D arrays with the same shape as z, or
    1D arrays with the
    appropriate dimensions.

    Keyword arguments:
    -----------------

    xr, xrange, yr, yrange
        plot ranges for x and y
    fill
        boolean parameter controlling whether to fill area between contours
        or not
    bar
        boolean parameter for plotting of the color-bar
    overplot
        boolean parameter for overplotting
    nlevels
        integer number of number of contours
    c_label
        boolean parameter for labeling or not-labeling each contour with
        its value
    """
    # Initialize x and y if these are not provided:
    if x is None or y is None:
        if z.ndim != 2:
            raise Exception("The 2D array is required")
        x, y = np.mgrid[0:z.shape[0], 0:z.shape[1]]
    else:
        if x.ndim == 1:
            x_new = x[:, np.newaxis] * (y * 0 + 1)
        else:
            x_new = x
        if y.ndim == 1:
            y_new = ((x * 0 + 1) * y[:, np.newaxis]).transpose()
        else:
            y_new = y
        x = x_new
        y = y_new

# Define position of this plot:
    if not noerase and not overplot:
        plt.gcf().clf()
        if position is not None:
            mypos = position[:]
            mypos[2] = position[2] - position[0]
            mypos[3] = position[3] - position[1]
            plt.axes(mypos)
    axis = plt.gca()

    # Rescaling of the axes:
    if xlog:
        axis.set_xscale('log')
    if ylog:
        axis.set_yscale('log')
    if zlog:
        z = np.log10(z)

# Setup axis ranges:
    if xr is not None:
        xrange = xr
    if yr is not None:
        yrange = yr
    if zr is not None:
        zrange = zr

    xrange = xrange or [x.min(), x.max()]
    yrange = yrange or [y.min(), y.max()]
    zrange = zrange or [z.min(), z.max()]

    # Setup levels for contour plot:
    if levels is None:
        zmin = zrange[0]
        zmax = zrange[1]
        levels = np.linspace(zmin, zmax, nlevels)

# Setup frame thickness:
#    axis.frame.set_linewidth(thick)

# Setup x-, y-titles and the main title:

    if xtitle is not None:
        axis.set_xlabel(xtitle)
    if ytitle is not None:
        axis.set_ylabel(ytitle)
    if title is not None:
        plt.title(title)

    # Assume autoscale is off:
    axis.set_autoscale_on(False)

    # Setup format of the major ticks:
    if xticklabel is not None:
        xFormatter = matplotlib.ticker.FormatStrFormatter(xticklabel)
        axis.xaxis.set_major_formatter(xFormatter)
    if yticklabel is not None:
        yFormatter = matplotlib.ticker.FormatStrFormatter(yticklabel)
        axis.yaxis.set_major_formatter(yFormatter)
    if xaxis_formatter is not None:
        axis.xaxis.set_major_formatter(xaxis_formatter)
    if yaxis_formatter is not None:
        axis.yaxis.set_major_formatter(yaxis_formatter)

    if not overplot:
        plt.gca().minorticks_on()

# Make a filled contour plot:
    if fill:
        cset1 = axis.contourf(x,
                              y,
                              z,
                              levels,
                              cmap=plt.cm.get_cmap(cm, nlevels))

# Add contour lines:
    if c_levels is None:
        c_levels = levels  # [0:len(levels):int(nlevels/12)]
    cset2 = axis.contour(x,
                         y,
                         z,
                         c_levels,
                         colors=c_color,
                         linewidths=c_thick,
                         hold='on')
    if label is not None:
        cset2.collections[0].set_label(label)
# Do not display dashed contours for negative values:
    for c in cset2.collections:
        c.set_linestyle(c_line)

# Add value labels on contour lines:
    if zticklabel is not None:
        zFormatter = matplotlib.ticker.FormatStrFormatter(zticklabel)
    else:
        zFormatter = None

    if c_label:
        if zticklabel is not None:
            args = {'fmt': zticklabel}
        else:
            args = {}
        axis.clabel(cset2, c_levels, inline=1, fontsize=c_charsize, **args)


# Do we need a color bar?:
    if fill & bar:
        #        matplotlib.rcParams['ytick.labelsize']=c_charsize
        if int(''.join((matplotlib.__version__).split('.')[:2])) >= 11:
            kw = {'use_gridspec': True}
        else:
            kw = {}
        plt.colorbar(cset1,
                     ticks=[np.min(levels), np.max(levels)],
                     fraction=bar_fraction,
                     format=zFormatter,
                     **kw)

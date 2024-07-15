from collections.abc import Iterable
import idlplot
import matplotlib.pyplot as plt
import numpy as np
"""
This module is a set wrappers around idlplot designed to make
the plots of subsets of the data: e.g.
plot(x,y,ind=ind)
instead of plot(x[ind],y[ind])
 
"""


def tvhist2d(a, b, *args, **kw):
    ind = kw.get('ind')

    if ind is None:
        return idlplot.tvhist2d(a, b, *args, **kw)
    else:
        weights = kw.get('weights')
        if weights is not None:
            kw['weights'] = kw['weights'][ind]
        del kw['ind']
        return idlplot.tvhist2d(a[ind], b[ind], *args, **kw)


def plothist(a, *args, **kw):
    ind = kw.get('ind')

    if ind is None:
        ret = idlplot.plothist(a, *args, **kw)
    else:
        weights = kw.get('weights')
        if weights is not None:
            if not np.isscalar(kw['weights']):
                kw['weights'] = kw['weights'][ind]
        del kw['ind']
        ret = idlplot.plothist(a[ind], *args, **kw)
    return ret


def plot(a, b=None, **kw):
    ind = kw.get('ind')

    if ind is None:
        idlplot.plot(a, b, **kw)
    else:
        del kw['ind']
        if b is not None:
            idlplot.plot(a[ind], b[ind], **kw)
        else:
            idlplot.plot(a[ind], None, **kw)


def plot_scatter(a, b, s=None, c=None, **kw):
    ind = kw.get('ind')

    if ind is None:
        idlplot.plot_scatter(a, b, s=s, c=c, **kw)
    else:
        del kw['ind']
        if c is not None and isinstance(c, Iterable):
            c = c[ind]
        if s is not None and isinstance(s, Iterable):
            s = s[ind]
        idlplot.plot_scatter(a[ind], b[ind], s=s, c=c, **kw)


def oplot(a, b=None, **kw):
    ind = kw.get('ind')

    if ind is None:
        idlplot.oplot(a, b, **kw)
    else:
        del kw['ind']
        if b is not None:
            idlplot.oplot(a[ind], b[ind], **kw)
        else:
            idlplot.oplot(a[ind], **kw)


def errorfixer(var, ind):
    var = np.asarray(var)
    if var.ndim == 2 and var.shape[0] == 2:
        var1 = [var[0][ind], var[1][ind]]
    else:
        var1 = var[ind]
    return var1


def ploterror(a, b, c, *args, **kw):
    ind = kw.get('ind')

    if ind is None:
        idlplot.ploterror(a, b, c, *args, **kw)
    else:
        del kw['ind']
        ll = len(args)
        args1 = [None] * ll
        c1 = errorfixer(c, ind)

        for i in range(ll):
            args1[i] = errorfixer(args[i], ind)
        idlplot.ploterror(a[ind], b[ind], c1, *args1, **kw)


def scatter(a, b, c=None, s=None, *args, **kw):
    ind = kw.get('ind')
    if ind is None:
        plt.scatter(a, b, c=c, s=s, *args, **kw)
    else:
        del kw['ind']
        if c is not None:
            c = c[ind]
        if s is not None:
            s = s[ind]
        plt.scatter(a[ind], b[ind], c=c, s=s, *args, **kw)


tvhist2d.__doc__ = idlplot.tvhist2d.__doc__
plot.__doc__ = idlplot.plot.__doc__
plot_scatter.__doc__ = idlplot.plot_scatter.__doc__
oplot.__doc__ = idlplot.oplot.__doc__
ploterror.__doc__ = idlplot.ploterror.__doc__
plothist.__doc__ = idlplot.plothist.__doc__

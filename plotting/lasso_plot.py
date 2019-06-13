# Copyright (C) 2010-2019 Sergey Koposov
# This file is part of astrolibpy
#
#    astrolibpy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#   astrolibpy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with astrolibpy.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
from matplotlib.widgets import Lasso
import matplotlib.path as mplpa
from matplotlib.pyplot import gca
import numpy as np
import quick_hist

class lasso_plot:
    """ The class is designed to select the datapoints by drawing the region
    around them.
    Example:
    plt.plot(xs, ys) # first plot the data
    las = lasso_plot.lasso_plot(xs,ys)
    Now click on the plot and do not release the mouse button till
    you draw your region
    After that the variable las.mask will contain the boolean mask of the
    points inside the region and las.verts will contain the vertices of the
    polygon you've just drawn
    The option bins is helpful when the dataset is very large. Then
    the data is binned first before checking whether it is inside the
    contour
    plt.plot(xs, ys)
    las = lasso_plot.lasso_plot(xs,ys,bins=200)
    """

    def __init__(self, xs, ys, bins=None):
        self.axes = gca()
        self.canvas = self.axes.figure.canvas
        self.xys = (np.asarray(xs), np.asarray(ys))
        self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)
        self.mask = None
        self.verts = None
        self.bins = bins

    def __getstate__(self):
        """ Custom pickle method to get rid of canvas/axes/lasso objects"""
        state = self.__dict__.copy()
        del state['canvas']
        del state['axes']
        del state['lasso']
        return state

    def callback(self, verts):
        self.verts = np.array(verts)
        mask = self.inside(self.xys[0], self.xys[1], self.bins)
        self.mask = mask
        self.canvas.draw_idle()
        self.canvas.widgetlock.release(self.lasso)
        self.canvas.mpl_disconnect(self.cid)
        del self.xys

    def inside(self, xs, ys, bins=None):
        """ Check if points xs,ys are inside the hand-selected mask """
        return inside(xs, ys, self.verts, bins)

    def onpress(self, event):
        if self.canvas.widgetlock.locked():
            return
        if event.inaxes is None:
            return
        self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata),
                           self.callback)
        # acquire a lock on the widget drawing
        self.canvas.widgetlock(self.lasso)


def inside(xs, ys, verts, bins=None):
    """ Check if points xs,ys are inside a mask specified by an array
    of vertices Nx2"""
    xys = np.array([xs, ys]).T

    path = mplpa.Path(verts)
    if bins is None:
        mask = path.contains_points(xys)
        ind = np.nonzero(mask)[0]
    else:
        if hasattr(bins, '__iter__'):
            if len(bins) == 2:
                nbins = bins
            else:
                raise Exception(
                    'The bins parameter must be a scalar or a tuple with 2 elements'
                )
        else:
            nbins = [bins, bins]
        minx, maxx = verts[:, 0].min(), verts[:, 0].max()
        miny, maxy = verts[:, 1].min(), verts[:, 1].max()
        hh, pos = quick_hist.quick_hist(xys.T,
                                        nbins=nbins,
                                        range=[[minx, maxx], [miny, maxy]],
                                        getPos=True)
        xbincens = np.linspace(minx, maxx, bins + 1, True)
        ybincens = np.linspace(miny, maxy, bins + 1, True)

        xbincens = .5 * (xbincens[1:] + xbincens[:-1])
        ybincens = .5 * (ybincens[1:] + ybincens[:-1])
        xbincens = xbincens[:, None] + ybincens[None, :] * 0
        ybincens = ybincens[None, :] + xbincens[:, None] * 0
        xybincens = np.array([xbincens.flatten(), ybincens.flatten()])

        mask = path.contains_points(xybincens.T)
        mask = mask[pos] & (pos >= 0)
    return mask

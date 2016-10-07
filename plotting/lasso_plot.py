# Copyright (C) 2010 Sergey Koposov
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
	"""
	def __init__(self, xs, ys, bins=None):
		self.axes = gca()
		self.canvas = self.axes.figure.canvas
		self.xys = np.array([xs,ys]).T#[d for d in zip(xs,ys)]
		fig = self.axes.figure
		self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)
		self.ind = None
		self.mask = None
		self.verts = None
		self.bins = bins

	def callback(self, verts):
		verts = np.array(verts)
		self.path = mplpa.Path(verts)
		if self.bins is None:
			mask = self.path.contains_points(self.xys)
			ind = np.nonzero(mask)[0]

		else:
			minx,maxx=verts[:,0].min(),verts[:,0].max()
			miny,maxy=verts[:,1].min(),verts[:,1].max()
			hh,pos= quick_hist.quick_hist(self.xys.T, nbins = [self.bins]*2,
					range=[[minx,maxx],[miny,maxy]],getPos=True)
			xbincens = np.linspace(minx,maxx,self.bins+1,True)
			ybincens = np.linspace(miny, maxy, self.bins+1,True)

			xbincens = .5*(xbincens[1:]+xbincens[:-1])
			ybincens = .5*(ybincens[1:]+ybincens[:-1])
			xbincens = xbincens[:,None] + ybincens[None,:]*0
			ybincens = ybincens[None,:] + xbincens[:,None]*0
			xybincens = np.array([xbincens.flatten(),ybincens.flatten()])	
			
			mask = self.path.contains_points(xybincens.T)
			mask = mask[pos] & (pos>=0)
			ind = mask	
		self.canvas.draw_idle()
		self.canvas.widgetlock.release(self.lasso)
		self.canvas.mpl_disconnect(self.cid)		
		del self.lasso
		del self.xys
		self.verts = verts
		self.ind = ind
		self.mask = mask

	def inside(self, xs,ys):
		if self.bins is None:
			tmpxys = zip(xs,ys)
			return self.path.contains_points(tmpxys)
		else:
			pass

	def onpress(self, event):
		if self.canvas.widgetlock.locked(): return
		if event.inaxes is None: return
		self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata), self.callback)
		# acquire a lock on the widget drawing
		self.canvas.widgetlock(self.lasso)
	

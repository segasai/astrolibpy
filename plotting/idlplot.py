# Copyright (C) 2009-2010 Sergey Koposov
# This file is part of astrolibpy
#
#	astrolibpy is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	astrolibpy is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with astrolibpy.  If not, see <http://www.gnu.org/licenses/>.


import matplotlib.pyplot as plt
import numpy
import scipy
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
import scipy.ndimage.filters
from matplotlib.pyplot import draw_if_interactive

import matplotlib
import types, sys
import warnings

# this module is by default in interactive regime 
plt.ion()

def listToArr(x):
	if isinstance(x, types.ListType):
		return numpy.array(x)
	else:
		return x

def listToArrFlat(x):
	if isinstance(x, types.ListType):
		return numpy.array(x).flatten()
	else:
		return x.flatten()

def get_marker(ps, linestyle):
	"""
	Wrapper for point markers which understand idl-like ps options
	(e.g. ps=3 for the point, ps=4 for the diamond)
	"""
	outlinestyle=' '
	markerHash={3: '.', 4:'D', 7:'+',2:'*',6:'s'}
	try:
		marker = markerHash[ps]
	except KeyError:
		if isinstance(ps,types.StringType):
			marker=ps
		else:
			marker=' '
			if linestyle is not None:
				outlinestyle=linestyle
			else:
				outlinestyle='-'
	return (marker, outlinestyle)

def exceptionDecorator(func):
	def wrapper(*args, **kwargs):
		try:
			isInteractive = plt.isinteractive()

			# switch to non-interactive mode
			matplotlib.interactive(False)

			ret = func(*args, **kwargs)

			matplotlib.interactive(isInteractive)

			draw_if_interactive()
			return ret
		except Exception, exc:
			# switch back
			matplotlib.interactive(isInteractive)

			einfo = sys.exc_info()
			raise einfo[0], einfo[1], einfo[2]

	return wrapper

def plothist(x, bin=None, nbins=None, xrange=None, yrange=None, min=None,
			max=None, overplot=False, color='black', xlog=False, ylog=False,
			nan=False, weights=None, norm=False, kernel=None, retpoints=False,
			adaptive=False, adaptive_thresh=30, adaptive_depth=[2,10], **kw):
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
		computed histogram.
		If yes the tuple with two arrays (bin centers, Number of points in bins) 
		is returned
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
	"""
	if nan:
		ind = numpy.isfinite(x)
		if weights is not None:
			ind = numpy.isfinite(weights)&ind
		dat = x[ind]
		if weights is not None:
			weights =weights[ind]
	else:
		dat = x
	if min is None:
		min = numpy.min(dat)
	if max is None:
		max = numpy.max(dat)
	
	if bin is None and nbins is None:
		nbins = 100
		bin = (max - min) * 1. / nbins
	elif nbins is None:
		nbins = int((max - min) * 1. / bin)
	elif bin is None:
		bin = (max - min) * 1. / nbins
	else:
		warnings.warn("both bin= and nbins= keywords were specified in the plothist call",RuntimeWarning)
		pass
		# if both nbins and bin are defined I don't do anything 
		# it may be non-intuitive if kernel option is used, because
		# it uses both nbins and bin options
	if kernel is None:
		if adaptive is None:
			hh, loc = numpy.histogram(dat, range=(min, max), bins=nbins, weights=weights)
		else:
			import adabinner
			hh, loc = adabinner.hist(dat, xmin=min, xmax=max, hi=adaptive_depth,
						thresh=adaptive_thresh)
		hh1 = numpy.zeros(2*len(hh)+2)
		loc1 = numpy.zeros_like(hh1)
		hh1[1:-1:2]=hh
		hh1[2::2]=hh
		loc1[1:-1:2]=loc[:-1]
		loc1[2::2]=loc[1:]
		loc1[-1]=loc1[-2]
		loc1[0]=loc[0]
	else:
		loc1=numpy.linspace(min,max,nbins*5)
		import statistics
		hh1=statistics.pdf( dat,loc1,h=bin/2.,kernel=kernel)*bin*len(dat)
	if overplot:
		func = oplot 
	else:
		func = plot
	if norm:
		hh1=hh1*1./hh1.max()
	func(loc1, hh1, ps=0, color=color, xrange=xrange, yrange=yrange,
		xlog=xlog, ylog=ylog, **kw)
	if retpoints:
		return 0.5*(loc[1:]+loc[:-1]),hh

@exceptionDecorator	 
def plot (arg1, arg2=None, xrange=None, yrange=None, ps=0, thick=1, xtitle=None, ytitle=None,
		color='black', noerase=False, overplot=False,position=None, ylog=False,
		xlog=False, xr=None, yr=None, title=None, label=None, nodata=False,
		linestyle=None, markersize=None, xaxis_formatter=None,
		yaxis_formatter=None, autoscalex=False, autoscaley=False,
		markerfacecolor=None,markeredgecolor=None):
	""" Plot your data in an IDL-way
		Example:
		plot(x,y,xrange=[0,39],yrange=[-1,10],ps=4,xtitle="X",\
			color='black',position=[0.1,0.1,0.9,0.9], xlog=True)
	"""

	if arg2 is None:
		y=listToArrFlat(arg1)
		x=numpy.arange(len(y))
	else:
		x=listToArrFlat(arg1)
		y=listToArrFlat(arg2)
		
	if not noerase:
		plt.gcf().clf()	
	if position is not None:
		mypos=position[:]
		mypos[2]=position[2]-position[0]
		mypos[3]=position[3]-position[1]
		plt.axes(mypos)
	axis = plt.gca()
	if xlog:
		axis.set_xscale('log',subx=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
	if ylog:
		axis.set_yscale('log',suby=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
	if xaxis_formatter is not None:
		axis.xaxis.set_major_formatter(xaxis_formatter)
	if yaxis_formatter is not None:
		axis.yaxis.set_major_formatter(yaxis_formatter)
	
	marker, linestyle = get_marker(ps, linestyle)
	if xr is not None:
		xrange=xr
	if yr is not None:
		yrange=yr

	if xrange is None and yrange is None:
		ind = numpy.isfinite(x)
		if not ind.any():
			xrange=[0,1]
		else:
			xrange=[numpy.min(x[ind]),numpy.max(x[ind])]
		ind = numpy.isfinite(y)
		if not ind.any():
			yrange=[0,1]
		else:
			yrange=[numpy.min(y[ind]),numpy.max(y[ind])]
		del ind
	elif xrange is None and yrange is not None:
		ind=(y<numpy.maximum(yr[1],yr[0])) & (y>numpy.minimum(yr[0],yr[1])) & numpy.isfinite(x)
		if ind.any():
			xrange=[numpy.min(x[ind]),numpy.max(x[ind])]
		else:
			xrange=[numpy.min(x),numpy.max(x)]
		del ind
	elif xrange is not None and yrange is None:
		ind=(x<numpy.maximum(xr[1],xr[0])) & (x>numpy.minimum(xr[0],xr[1])) & numpy.isfinite(y)
		if ind.any():
			yrange=[numpy.min(y[ind]),numpy.max(y[ind])]
		else:
			yrange=[numpy.min(y),numpy.max(y)]
		del ind
	if len(yrange)!=2 or len(xrange)!=2:
		raise ValueError("Wrong xrange or yrange")
	if not overplot:
		axis.minorticks_on()
	if xtitle is not None:
		axis.set_xlabel(xtitle)
	if ytitle is not None:
		axis.set_ylabel(ytitle)
		
	axis.set_autoscalex_on(autoscalex)
	axis.set_autoscaley_on(autoscaley)
	if not overplot:
		axis.axis(numpy.concatenate((xrange,yrange)))
	if title is not None:
		plt.title(title)
	if not nodata:
		if markersize is None:
			axis.plot(x, y, marker=marker, linestyle=linestyle,
							linewidth=thick, color=color, label=label,
							markerfacecolor=markerfacecolor,
							markeredgecolor=markeredgecolor)
		else:
			axis.plot(x, y, marker=marker, linestyle=linestyle,
							linewidth=thick, color=color, label=label,
							markersize=markersize,
							markerfacecolor=markerfacecolor,
							markeredgecolor=markeredgecolor)	
	
def oplot (x, y=None, **kw):
	"""
	Overplot your data
	
	Example:
	>> oplot(x,2+y/10.,ps=3,color='blue')
	"""

	plot (x,y, noerase=True, overplot=True, **kw)

@exceptionDecorator
def ploterror (x, y, err0, err1=None, color='black', ps=0, ecolor='black',
				overplot=False, noerase=False, elinewidth=None, capsize=None,
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
		noerase=True
	if err1 is None:
		erry = listToArr(err0)
	else:
		erry = listToArr(err1)
		errx = listToArr(err0)
	
	if kw.get('yr') is None:
		kw['yr'] = [(y-erry).min(),(y+erry).max()]
	plot (x, y, color=color, ps=ps, overplot=overplot, noerase=noerase, **kw)
	(marker, outlinestyle) = get_marker(ps, None)
	kw1 = {'ecolor':ecolor, 'marker':marker, 'color':color, 'linestyle':outlinestyle,
			'elinewidth':elinewidth} 
	if capsize is not None:
		kw1['capsize']=capsize
	if err1 is None:
		plt.gca().errorbar(x, y, erry, **kw1)
	else:
		plt.gca().errorbar(x, y, xerr=errx, **kw1)
		plt.gca().errorbar(x, y, yerr=erry, **kw1)

@exceptionDecorator
def tvaxis (image, xmin=None, xmax=None, ymin=None,ymax=None, xtitle="", ytitle="", title="",
			vmin=None, vmax=None, aspect="auto", xlog=False ,ylog=False,
			position=None, noerase=False, bar=False, bar_label='',
			bar_fraction=0.05, zlog=False, smooth=None, **kw):
	"""
	Display the 2D image with proper axes (similar to plt.imshow)
	Example:
	>> tvaxis(im,-20,10,-40,50)
	
	Keyword parameters:
	------------------
	xmin,xmax,ymin,ymax
		the ranges for x,y where the histogram is constructed.
		These params can be specified not only as keywords but also as normal arguments
	vmin,vmax
		the ranges of intensities shown on the plot
	smooth
		if not None this parameter controls additional smoothing of the 2D histogram
	bar
		boolean parameter for switching on/off the plotting of the color-bar
	
	"""

	if xlog:
		plt.gca().set_xscale('log')
	if ylog:
		plt.gca().set_yscale('log')
	if not noerase:
		plt.clf()

	if position is not None:
		mypos=position[:]
		mypos[2]=position[2]-position[0]
		mypos[3]=position[3]-position[1]
		plt.axes(mypos)
	if xmin is None:
		xmin = 0
	if ymin is None:
		ymin = 0
	if xmax is None:
		xmax = image.shape[0]
	if ymax is None:
		ymax = image.shape[1]
	
	im = image.T
	if smooth is not None:
		im = scipy.ndimage.filters.gaussian_filter(im, [smooth,smooth])
	if zlog:
		norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
	else:
		norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)		
	
	axim = plt.imshow(im, extent=(xmin, xmax, ymin, ymax), vmin=vmin, vmax=vmax, 
					aspect=aspect, norm=norm, **kw)

	plt.gca().set_xlabel(xtitle)
	plt.gca().set_ylabel(ytitle)
	plt.gca().minorticks_on()

	if title is not None:
		plt.title(title)
	if bar:
		cb=plt.colorbar(fraction=bar_fraction)
		cb.set_label(bar_label)

	return axim

@exceptionDecorator
def tvhist2d (x,y, xmin=None, xmax=None, ymin=None, ymax=None,
				vmin=None, vmax=None, bins=[100,100], xtitle="",
				ytitle="", noerase=False, weights=None, zlog=False,
				xflip=False, yflip=False, bar=False, bar_label='',
				bar_fraction=0.05, smooth=None, quick=False,
				cmap='gray_r', normx=False, normy=False,
				xlog=False, ylog=False, **kw):
	""" Plot the 2D histogram of the data
	Example:
	>> tvhist2d(xs,ys,bins=[30,30])
	
	>> tvhist2d(xs,ys,0,10,-1,2,bins=[30,30])
	
	Keyword arguments:
	-----------------
	xmin,xmax,ymin,ymax
		the ranges for x,y where the histogram is constructed.
		These params can be specified not only as keywords but also as normal arguments
		>> tvhist2d(xs,ys,xmin=0,ymin=10,ymin=-1,ymax=2,bins=[30,30])
		>> tvhist2d(xs,ys,0,10,-1,2,bins=[30,30])
	vmin,vmax
		the ranges of intensities shown on the plot
	bins
		the list of two integers specifying how many bins in x,y you want
	smooth
		if not None this parameter controls additional smoothing of the 2D histogram
	bar
		boolean parameter for switching on/off the plotting of the color-bar
	xflip, yflip
		boolean parameters allowing to flip x,y axes
	normx, normy
		boolean params controlling the normalization of the histogram along X or Y axes
		in such way that the brightest pixel value will be the same for each row/column
	"""

	x1 = listToArrFlat(x)
	y1 = listToArrFlat(y)
	ind = numpy.isfinite(x1) & numpy.isfinite(y1)

	if xmin is None:
		xmin = x1[ind].min()
	if ymin is None:
		ymin = y1[ind].min()
	if xmax is None:
		xmax = x1[ind].max()
	if ymax is None:
		ymax = y1[ind].max()

	range1 = (xmin, xmax, ymin, ymax)
	if xlog is True:
		x1=numpy.log10(x1)
		xmin,xmax=numpy.log10(xmin),numpy.log10(xmax)
	if ylog is True:
		y1=numpy.log10(y1)
		ymin,ymax=numpy.log10(ymin),numpy.log10(ymax)
	range = [[ymin, ymax],[xmin, xmax]]
	ind = numpy.isfinite(x1) & numpy.isfinite(y1)

	if not quick:
		hh, yedges, xedges = scipy.histogram2d(y1[ind], x1[ind], range=range,
												bins=bins, weights=weights)
	else:
		import quick_hist
		hh = quick_hist.quick_hist((y1[ind], x1[ind]), range=range, nbins=bins,
								weights=weights)
	if normx:
		hh = hh*1./numpy.maximum(hh.sum(axis=0),1)[numpy.newaxis,:]
	if normy:
		hh = hh*1./numpy.maximum(hh.sum(axis=1),1)[:,numpy.newaxis]

	if xflip:
		range1 = (range1[1], range1[0], range1[2], range1[3])
		hh = numpy.fliplr(hh)
	if yflip:
		range1 = (range1[0], range1[1], range1[3], range1[2])
		hh = numpy.flipud(hh)

	if not noerase:
		plt.gcf().clf()

	plt.gca().set_xlabel(xtitle)
	plt.gca().set_ylabel(ytitle)
	if smooth is not None:
		hh = scipy.ndimage.filters.gaussian_filter(hh, [smooth, smooth])
	if zlog:
		norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
	else:
		norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)		
	axim=plt.imshow(hh, extent=range1, aspect='auto', interpolation='nearest',
					cmap=cmap, norm=norm, **kw)
	if bar:
		cb=plt.colorbar(fraction=bar_fraction)
		cb.set_label(bar_label)
	if xlog:
		plt.gca().set_xscale('log')
	if ylog:
		plt.gca().set_yscale('log')

	return axim

@exceptionDecorator
def contour (z, x=None, y=None, xrange=None, yrange=None, zrange=None,
		xr=None, yr=None, zr=None, title="", xtitle="", ytitle="",
		position=None, xlog=False, ylog=False, zlog=False, xticklabel = None,
		yticklabel = None, zticklabel = None, levels=None, nlevels=256,
		c_color="black", cm ="jet", c_line="solid", c_levels=None,
		c_charsize=12.0, c_thick=1, thick=1, font="monospace",
		weight="normal", charsize=14.0, bar=True, fill=True, overplot=False,
		noerase=False, c_label=False, bar_fraction=0.05, xaxis_formatter=None,
		yaxis_formatter=None):
	"""
	Plot the contours of the 2d array. 
	Example:
	>> contour(z)
	if you have the x and y coordinates of your array then you can use them
	>> contour(z, x, y)
	In that case  x,y can be either 2D arrays with the same shape as z, or 1D arrays with the
	appropriate dimensions.
	
	Keyword arguments:
	-----------------
	
	xr, xrange, yr, yrange
		plot ranges for x and y
	fill
		boolean parameter controlling whether to fill area between contours or not
	bar
		boolean parameter for plotting of the color-bar
	overplot
		boolean parameter for overplotting
	nlevels
		integer number of number of contours
	c_label
		boolean parameter for labeling or not-labeling each contour with its value
	"""
	# Initialize x and y if these are not provided:
	if x is None or y is None:
		if z.ndim!=2:
			raise Exception("The 2D array is required")
		x, y = numpy.mgrid[0:z.shape[0],0:z.shape[1]]
	else:
		if x.ndim==1:
			x_new=x[:,numpy.newaxis]*(y*0+1)
		else:
			x_new=x
		if y.ndim==1:
			y_new=((x*0+1)*y[:,numpy.newaxis]).transpose()
		else:
			y_new=y
		x=x_new
		y=y_new

# Define position of this plot:
	if not noerase and not overplot:
		plt.gcf().clf()	
		if position is not None:
			mypos=position[:]
			mypos[2]=position[2]-position[0]
			mypos[3]=position[3]-position[1]
			plt.axes(mypos)
	axis = plt.gca()

# Rescaling of the axes:		
	if xlog:
	 	axis.set_xscale('log')
	if ylog:
		axis.set_yscale('log')
	if zlog:
		z = numpy.log10(z)

# Setup axis ranges:			
	if xr is not None:
		xrange=xr
	if yr is not None:
		yrange=yr
	if zr is not None:
		zrange=zr

	xrange = xrange or [x.min(),x.max()]
	yrange = yrange or [y.min(),y.max()]
	zrange = zrange or [z.min(),z.max()]		
		
# Setup levels for contour plot:		
	if levels is None:
		zmin = zrange[0]
		zmax = zrange[1]
		levels = numpy.linspace(zmin,zmax,nlevels)

# Setup frame thickness:
#	axis.frame.set_linewidth(thick) 

# Setup x-, y-titles and the main title:
	axis.set_xlabel(xtitle)
	axis.set_ylabel(ytitle)
	plt.title(title)

# Assume autoscale is off:	
	axis.set_autoscale_on(False)
	
# For a new plot, we have to initialize axes:	
	if not overplot:
		axis.axis(xrange+yrange) 

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
		cset1=axis.contourf(x, y, z, levels, cmap=plt.cm.get_cmap(cm, nlevels))

# Add contour lines:
	if c_levels is None:
		c_levels = levels[0:len(levels):int(nlevels/12)]
	cset2 = axis.contour(x, y, z, c_levels, colors = c_color,linewidths=c_thick,hold='on')

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
			args = {'fmt':zticklabel}
		else:
			args = {}
		cset3 = axis.clabel(cset2, c_levels, inline=1, fontsize=c_charsize, **args)
				
# Do we need a color bar?:
	if fill & bar:
#		matplotlib.rcParams['ytick.labelsize']=c_charsize				
		plt.colorbar(cset1, ticks=[numpy.min(levels), numpy.max(levels)],#, shrink = 0.87, aspect=15, 
			fraction=bar_fraction, format=zFormatter)

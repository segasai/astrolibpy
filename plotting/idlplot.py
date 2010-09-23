# Copyright (C) 2009-2010 Sergey Koposov
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


import matplotlib.pyplot as plt
import numpy
import scipy
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
import scipy.ndimage.filters

import matplotlib
import types

plt.ion()

def get_marker(ps, linestyle):
	outlinestyle=' '
	if ps==3:
		marker='.'
	elif ps==4: 
		marker='D'
	elif ps==7: 
		marker='+'
	elif ps==2: 
		marker='*'
	elif ps==6: 
		marker='s'
	else:
		if isinstance(ps,types.StringType):
			marker=ps
		else:
			marker=' '
			if linestyle is not None:
				outlinestyle=linestyle
			else:
				outlinestyle='-'
	return (marker, outlinestyle)

def plothist(x,bin=None, xrange=None, yrange=None, min=None, max=None,
			overplot=False,color='black', xlog=False, ylog=False,
			nan=False, weights=None, norm=False, kernel=None, **kw):
	if nan:
		ind = numpy.isfinite(x)
		dat = x[ind]
	else:
		dat = x
	if min is None:
		min=numpy.min(dat)
	if max is None:
		max=numpy.max(dat)
	if bin is None:
		nbin=100
		bin = (max-min)*1./nbin
	else:
		nbin=(max-min)/bin
	if kernel is None:
		hh, loc = numpy.histogram(dat, range=(min, max), bins=nbin, weights=weights)
		hh1 = numpy.zeros(2*len(hh)+2)
		loc1 = numpy.zeros_like(hh1)
		hh1[1:-1:2]=hh
		hh1[2::2]=hh
		loc1[1:-1:2]=loc[:-1]
		loc1[2::2]=loc[1:]
		loc1[-1]=loc1[-2]
		loc1[0]=loc[0]
	else:
		loc1=numpy.linspace(min,max,nbin*5)
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
	 
def plot (arg1, arg2=None, xrange=None, yrange=None, ps=0, thick=1, xtitle=None, ytitle=None,
		color='black', noerase=False, overplot=False,position=None, ylog=False,
		xlog=False, xr=None, yr=None, title=None, label=None, nodata=False,
		linestyle=None, markersize=None, xaxis_formatter=None,
		yaxis_formatter=None, autoscalex=False, autoscaley=False,
		markerfacecolor=None):
	""" Plot your data in an IDL-way
		Example:
		plot(x,y,xrange=[0,39],yrange=[-1,10],ps=4,xtitle="X",\
			color='black',position=[0.1,0.1,0.9,0.9], xlog=True)
	"""
	listtoarr = lambda x: numpy.array(x) if isinstance(x, types.ListType) else x

	if arg2 is None:
		y=listtoarr(arg1)
		x=numpy.arange(len(y))
	else:
		x=listtoarr(arg1)
		y=listtoarr(arg2)
	if x.ndim !=1:
		x=x.flatten()
	if y.ndim !=1:
		y=y.flatten()
		
	if not noerase:
		plt.gcf().clf()	
	if position is not None:
		mypos=position[:]
		mypos[2]=position[2]-position[0]
		mypos[3]=position[3]-position[1]
		plt.axes(mypos)
	
	if xlog:
		plt.gca().set_xscale('log',subx=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
	if ylog:
		plt.gca().set_yscale('log',suby=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
	if xaxis_formatter is not None:
		plt.gca().xaxis.set_major_formatter(xaxis_formatter)
	if yaxis_formatter is not None:
		plt.gca().yaxis.set_major_formatter(yaxis_formatter)
	
	marker, linestyle = get_marker(ps, linestyle)
	if xr is not None:
		xrange=xr
	if yr is not None:
		yrange=yr

	if xrange is None and yrange is None:
		ind = numpy.isfinite(x)
		xrange=[numpy.min(x[ind]),numpy.max(x[ind])]
		ind = numpy.isfinite(y)
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
		if not xlog:
			xminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])
			plt.gca().xaxis.set_minor_locator(xminorLocator)
		if not ylog:
			yminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])
			plt.gca().yaxis.set_minor_locator(yminorLocator)		

	if xtitle is not None:
		plt.gca().set_xlabel(xtitle)
	if ytitle is not None:
		plt.gca().set_ylabel(ytitle)
        
	plt.gca().set_autoscalex_on(autoscalex)
	plt.gca().set_autoscaley_on(autoscaley)
	if not overplot:
		plt.gca().axis(numpy.concatenate((xrange,yrange)))
	if title is not None:
		plt.title(title)
	if not nodata:
		if markersize is None:
			plt.gca().plot(x, y, marker=marker, linestyle=linestyle,
							linewidth=thick, color=color, label=label,
							markerfacecolor=markerfacecolor)
		else:
			plt.gca().plot(x, y, marker=marker, linestyle=linestyle,
							linewidth=thick, color=color, label=label,
							markersize=markersize,
							markerfacecolor=markerfacecolor)	
	if plt.isinteractive():
		plt.draw()
	
def oplot (x,y=None, **kw):
	""" Overplot your data in an IDL-way
		Example:
		idlplot.oplot(x,2+y/10.,ps=3,color='blue')
	"""

	plot (x,y, noerase=True, overplot=True, **kw)

def ploterror (x,y, err0, err1=None, color='black', ps=0, ecolor='black', overplot=False, 
				noerase=False, elinewidth=None, **kw):
	if overplot:
		noerase=True
	if err1 is None:
		erry = err0
	else:
		erry = err1
	if kw.get('yr') is None:
		kw['yr'] = [(y-erry).min(),(y+erry).max()]
	plot (x,y,color=color, ps=ps, overplot=overplot, noerase=noerase, **kw)
	(marker,outlinestyle)=get_marker(ps, None)
	if err1 is None:
		plt.gca().errorbar(x,y,err0,color=color,ecolor=ecolor,marker=marker,
						linestyle=outlinestyle,elinewidth=elinewidth)
	else:
		plt.gca().errorbar(x,y,xerr=err0,color=color,ecolor=ecolor,marker=marker,
						linestyle=outlinestyle,elinewidth=elinewidth)
		plt.gca().errorbar(x,y,yerr=err1,color=color,ecolor=ecolor,marker=marker,
						linestyle=outlinestyle,elinewidth=elinewidth)
	
	if plt.isinteractive():
		plt.draw()


def tvaxis (image, xmin=None, xmax=None, ymin=None,ymax=None, xtitle="", ytitle="", title="",
			vmin=None, vmax=None, aspect="auto", xlog=False ,ylog=False,
			position=None, noerase=False, bar=False, bar_label='',
			bar_fraction=0.05, zlog=False, smooth=None, **kw):

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
		xmin=0
	if ymin is None:
		ymin=0
	if xmax is None:
		xmax=image.shape[0]
	if ymax is None:
		ymax=image.shape[1]
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
	xminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])
	yminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])
#	seg=  mlines.Line2D(x, y, color=_color, linestyle=linestyle, 
#							marker=marker, axes=plt.gca() )

	plt.gca().xaxis.set_minor_locator(xminorLocator)
	plt.gca().yaxis.set_minor_locator(yminorLocator)		

	if title is not None:
		plt.title(title)
	if bar:
		cb=plt.colorbar(fraction=bar_fraction)
		cb.set_label(bar_label)
	if plt.isinteractive():
		plt.draw()
	return axim

def tvhist2d (x,y, xmin=None, xmax=None, ymin=None, ymax=None,
				vmin=None, vmax=None, bins=[100,100], xtitle="",
				ytitle="", noerase=False, weights=None, zlog=False,
				xflip=False, yflip=False, bar=False, bar_label='',
				bar_fraction=0.05, smooth=None, quick=False,
				cmap='gray_r', **kw):
	""" Plots the 2D histogram of the data"""
	if not noerase:
		plt.gcf().clf()
	x1 = x.flat
	y1 = y.flat
	ind = numpy.isfinite(x1) & numpy.isfinite(y1)
	if xmin is None:
		xmin = x1[ind].min()
	if ymin is None:
		ymin = y1[ind].min()
	if xmax is None:
		xmax = x1[ind].max()
	if ymax is None:
		ymax = y1[ind].max()
	range = [[ymin, ymax],[xmin, xmax]]
	range1 = (xmin, xmax, ymin, ymax)
	if not quick:
		hh, yedges, xedges = scipy.histogram2d(y1[ind], x1[ind], range=range,
		                                        bins=bins, weights=weights)
	else:
		import quick_hist
		hh = quick_hist.quick_hist((y1[ind], x1[ind]), range=range, nbins=bins,
                                weights=weights)
	if xflip:
		range1 = (range1[1], range1[0], range1[2], range1[3])
		hh = numpy.fliplr(hh)
	if yflip:
		range1 = (range1[0], range1[1], range1[3], range1[2])
		hh = numpy.flipud(hh)

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
	if plt.isinteractive():
		plt.draw()
	return axim


def contour (z, x=None, y=None, xrange=None, yrange=None, zrange=None,
		xr=None, yr=None, zr=None, title="", xtitle="", ytitle="",
		position=None, xlog=False, ylog=False, zlog=False, xticklabel = None,
		yticklabel = None, zticklabel = None, levels=None, nlevels=256,
		c_color="black", cm ="jet", c_line="solid", c_levels=None,
		c_charsize=12.0, c_thick=1, thick=1, font="monospace",
		weight="normal", charsize=14.0, bar=True, fill=True, overplot=False,
		noerase=False, c_label=True, bar_fraction=0.05, xaxis_formatter=None,
		yaxis_formatter=None):

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
	if not noerase:
		plt.gcf().clf()	
		if position is not None:
			mypos=position[:]
			mypos[2]=position[2]-position[0]
			mypos[3]=position[3]-position[1]
			plt.axes(mypos)

# Rescaling of the axes:		
	if xlog:
	 	plt.gca().set_xscale('log')
	if ylog:
		plt.gca().set_yscale('log')
	if zlog:
		z = numpy.log10(z)

# Setup axis ranges:			
	if xr is not None:
		xrange=xr
	if yr is not None:
		yrange=yr
	if zr is not None:
		zrange=zr

	if xrange is None:
		xrange=[x.min(),x.max()]
	if yrange is None:
		yrange=[y.min(),y.max()]
	if zrange is None:
		zrange=[z.min(),z.max()]		
		
# Setup levels for contour plot:		
	if levels is None:
		zmin = zrange[0]
		zmax = zrange[1]
		dz = zmax-zmin
		levels = zmin + dz/nlevels*numpy.arange(nlevels+1)	

# Setup frame thickness:
#	plt.gca().frame.set_linewidth(thick) 

# Setup x-, y-titles and the main title:
	plt.gca().set_xlabel(xtitle)
	plt.gca().set_ylabel(ytitle)
	plt.title(title)

# Assume autoscale is off:	
	plt.gca().set_autoscale_on(False)
	
# For a new plot, we have to initialize axes:	
	if not overplot:
		plt.gca().axis(xrange+yrange) 

# Setup format of the major ticks:
	if xticklabel is not None:
		xFormatter = matplotlib.ticker.FormatStrFormatter(xticklabel)
		plt.gca().xaxis.set_major_formatter(xFormatter)
	if yticklabel is not None:
		yFormatter = matplotlib.ticker.FormatStrFormatter(yticklabel)
		plt.gca().yaxis.set_major_formatter(yFormatter)
	if xaxis_formatter is not None:
		plt.gca().xaxis.set_major_formatter(xaxis_formatter)
	if yaxis_formatter is not None:
		plt.gca().yaxis.set_major_formatter(yaxis_formatter)

# Setup minor tickmarks:		
	xminorLocator = matplotlib.ticker.MaxNLocator(nbins=90, steps=[1, 2, 5, 10])
	yminorLocator = matplotlib.ticker.MaxNLocator(nbins=90, steps=[1, 2, 5, 10])	

	plt.gca().xaxis.set_minor_locator(xminorLocator)
	plt.gca().yaxis.set_minor_locator(yminorLocator)		
		
# Make a filled contour plot:		
	if fill:
		cset1=plt.gca().contourf(x, y, z, levels, cmap=plt.cm.get_cmap(cm, nlevels))

# Add contour lines:
	if c_levels is None:
		c_levels = levels[0:len(levels):int(nlevels/12)]
	cset2 = plt.gca().contour(x, y, z, c_levels, colors = c_color,linewidths=c_thick,hold='on')

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
		cset3 = plt.gca().clabel(cset2, c_levels, inline=1, fontsize=c_charsize, **args)
				
# Do we need a color bar?:
	if fill & bar:
#		matplotlib.rcParams['ytick.labelsize']=c_charsize				
		plt.colorbar(cset1, ticks=[numpy.min(levels), numpy.max(levels)],#, shrink = 0.87, aspect=15, 
			fraction=bar_fraction, format=zFormatter)
		
	if plt.isinteractive():
		plt.draw()
		


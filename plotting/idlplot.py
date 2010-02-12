import matplotlib.pyplot as plt
import numpy
import scipy
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
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
			nan=False, weights=None, norm=False, **kw):
	if min==None:
		min=numpy.nanmin(x)
	if max==None:
		max=numpy.nanmax(x)
	if bin==None:
		bin=100
	else:
		bin=(max-min)/bin
	if nan:
		ind = numpy.isfinite(x)
		dat = x[ind]
	else:
		dat = x
	hh, loc = numpy.histogram(dat, range=(min, max), bins=bin, weights=weights)
	hh1 = numpy.zeros(2*len(hh)+2)
	loc1 = numpy.zeros_like(hh1)
	hh1[1:-1:2]=hh
	hh1[2::2]=hh
	loc1[1:-1:2]=loc[:-1]
	loc1[2::2]=loc[1:]
	loc1[-1]=loc1[-2]
	loc1[0]=loc[0]
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
		yaxis_formatter=None):
	""" Plot your data in an IDL-way
		Example:
		plot(x,y,xrange=[0,39],yrange=[-1,10],ps=4,xtitle="X",\
			color='black',position=[0.1,0.1,0.9,0.9], xlog=True)
	"""
	listtoarr = lambda x: numpy.array(x) if isinstance(x, types.ListType) else x

	if arg2==None:
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
	if position!=None:
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
		xrange=[min(x[ind]),max(x[ind])]
		ind = numpy.isfinite(y)
		yrange=[min(y[ind]),max(y[ind])]
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

	if xtitle!=None:
		plt.gca().set_xlabel(xtitle)
	if ytitle!=None:
		plt.gca().set_ylabel(ytitle)
	plt.gca().set_autoscale_on(False)
	if not overplot:
		plt.gca().axis(numpy.concatenate((xrange,yrange)))
	if title!= None:
		plt.title(title)
	if not nodata:
		if markersize is None:
			plt.gca().plot(x, y, marker=marker, linestyle=linestyle,
							linewidth=thick, color=color, label=label)
		else:
			plt.gca().plot(x, y, marker=marker, linestyle=linestyle,
							linewidth=thick, color=color, label=label,
							markersize=markersize)	
	if plt.isinteractive():
		plt.draw()
	
def oplot (x,y=None, **kw):
	""" Overplot your data in an IDL-way
		Example:
		idlplot.oplot(x,2+y/10.,ps=3,color='blue')
	"""

	plot (x,y, noerase=True, overplot=True, **kw)

def ploterror (x,y, err, color='black', ps=0, ecolor='black', overplot=False, 
				noerase=False, **kw):
	if overplot:
		noerase=True
	if kw.get('yr') == None:
		kw['yr'] = [(y-err).min(),(y+err).max()]
	plot (x,y,color=color, ps=ps, overplot=overplot, noerase=noerase, **kw)
	(marker,outlinestyle)=get_marker(ps, None)
	plt.gca().errorbar(x,y,err,color=color,ecolor=ecolor,marker=marker,
						linestyle=outlinestyle)
	if plt.isinteractive():
		plt.draw()


def tvaxis (image, xmin=None, xmax=None, ymin=None,ymax=None, xtitle="", ytitle="", title="",
			vmin=None, vmax=None, aspect="auto", xlog=False ,ylog=False,
			position=None, noerase=False, bar=False, **kw):

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

	plt.imshow(image.transpose(), extent=(xmin, xmax, ymin, ymax), vmin=vmin, vmax=vmax, 
					aspect=aspect, **kw)

	plt.gca().set_xlabel(xtitle)
	plt.gca().set_ylabel(ytitle)
	xminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])
	yminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])
#	seg=  mlines.Line2D(x, y, color=_color, linestyle=linestyle, 
#							marker=marker, axes=plt.gca() )

	plt.gca().xaxis.set_minor_locator(xminorLocator)
	plt.gca().yaxis.set_minor_locator(yminorLocator)		

	if title!= None:
		plt.title(title)
	if bar:
		plt.colorbar()
	if plt.isinteractive():
		plt.draw()

def tvhist2d (x,y, xmin=None, xmax=None, ymin=None, ymax=None,
				bins=[100,100], xtitle="",
				ytitle="", noerase=False, weights=None,
				**kw):
	""" Plots the 2D histogram of the data"""
	if not noerase:
		plt.gcf().clf()
	if xmin is None:
		xmin = x.min()
	if ymin is None:
		ymin = y.min()
	if xmax is None:
		xmax = x.max()
	if ymax is None:
		ymax = y.max()
	range=[[ymin,ymax],[xmin,xmax]]
	range1=(xmin,xmax,ymin,ymax)

	hh,xedges,yedges=scipy.histogram2d(y,x,range=range, bins=bins, weights=weights)
	if range1 is None:
		range1=(yedges[0],yedges[-1],xedges[0],xedges[-1])
		print range1
	plt.gca().set_xlabel(xtitle)
	plt.gca().set_ylabel(ytitle)

	plt.imshow(-hh,extent=range1, aspect='auto', interpolation='nearest', **kw)
	if plt.isinteractive():
		plt.draw()



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
	if xr!=None:
		xrange=xr
	if yr!=None:
		yrange=yr
	if zr!=None:
		zrange=zr

	if xrange==None:
		xrange=[x.min(),x.max()]
	if yrange==None:
		yrange=[y.min(),y.max()]
	if zrange==None:
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
		


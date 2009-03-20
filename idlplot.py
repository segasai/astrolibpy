import matplotlib.pyplot as plt
import numpy
import scipy
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator

plt.ion()

#idlplot.plot(x,y,xrange=[0,39],yrange=[-1,10],ps=4,xtitle="X",\
#    color='black',position=[0.1,0.1,0.9,0.9], xlog=True)
#idlplot.oplot(x,2+y/10.,ps=3,xtitle="SX",color='blue')

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
		marker=' '
		outlinestyle='-'
	return (marker, outlinestyle)

def plothist(x,bin=None, xrange=None, yrange=None, min=None, max=None,
			overplot=False,color='black', xlog=False, ylog=False, **kw):
	if min==None:
		min=x.min()
	if max==None:
		max=x.max()
	if bin==None:
		bin=100
	else:
		bin=(max-min)/bin
	hh = numpy.histogram(x,range=(min,max),bins=bin)
	loc=hh[1]
	hh=hh[0]
	hh1 = numpy.zeros(2*len(hh)+2)
	loc1 = numpy.zeros(2*len(hh)+2)
	for a in range(len(hh)):
		hh1[2*a+1]=hh[a]
		hh1[2*a+2]=hh[a]
		loc1[2*a+1]=loc[a]
		loc1[2*a+2]=loc[a+1]
	loc1[-1]=loc1[-2]
	loc1[0]=loc[0]
	if overplot:
		tmp=oplot 
	else:
		tmp=plot
	tmp(loc1,hh1,ps=0,color=color,xrange=xrange,yrange=yrange,
		xlog=xlog,ylog=ylog, **kw)
    
def plot (arg1, arg2=None, xrange=None, yrange=None, ps=0, thick=1, xtitle=None, ytitle=None,
		color='black', noerase=False, overplot=False,position=None, ylog=False,
		xlog=False, xr=None, yr=None, title=None):
	""" Plot your data in an IDL-way
		Example:
		plot(x,y,xrange=[0,39],yrange=[-1,10],ps=4,xtitle="X",\
			color='black',position=[0.1,0.1,0.9,0.9], xlog=True)
	"""
	if arg2==None:
		y=arg1
		x=numpy.arange(len(y))
	else:
		x=arg1
		y=arg2
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
		plt.gca().set_xscale('log')
	if ylog:
		plt.gca().set_yscale('log')
	( marker, linestyle) = get_marker(ps, None)		                    
	if xr!=None:
		xrange=xr
	if yr!=None:
		yrange=yr

	if xrange is None and yrange is None:
		xrange=[min(x),max(x)]
		yrange=[min(y),max(y)]
	elif xrange is None and yrange is not None:
		ind=numpy.where((y<yr[1]) & (y>yr[0]))[0]
		if len(ind)!=0:
			xrange=[numpy.min(x[ind]),numpy.max(x[ind])]
		else:
			xrange=[numpy.min(x),numpy.max(x)]		
	elif xrange is not None and yrange is None:
		ind=numpy.where((x<xr[1]) & (x>xr[0]))[0]
		if len(ind)!=0:
			yrange=[numpy.min(y[ind]),numpy.max(y[ind])]
		else:
			yrange=[numpy.min(y),numpy.max(y)]		
		
	xminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])
	yminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])
#	seg=  mlines.Line2D(x, y, color=_color, linestyle=linestyle, 
#							marker=marker, axes=plt.gca() )

	plt.gca().xaxis.set_minor_locator(xminorLocator)
	plt.gca().yaxis.set_minor_locator(yminorLocator)		
	if xtitle!=None:
		plt.gca().set_xlabel(xtitle)
	if ytitle!=None:
		plt.gca().set_ylabel(ytitle)
	plt.gca().set_autoscale_on(False)
	if not overplot:
		plt.gca().axis(xrange+yrange)#,xautcoscale_on=False)
	if title!= None:
		plt.title(title)

	plt.gca().plot(x,y,marker=marker,linestyle=linestyle,linewidth=thick,color=color)
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
	plt.draw()		


def tvaxis (image, xmin, xmax,ymin,ymax):
	pass

def tvhist2d (x,y, xmin=None, xmax=None, ymin=None, ymax=None,
				bins=[100,100], xtitle="",
				ytitle="", noerase=False, **kw):
	""" Plots the 2D histogram of the data"""
	if not noerase:
		plt.gcf().clf()
	if xmin==None or xmax==None or ymin==None or ymax==None:
		range=None
		range1=None
	else:
		range=[[ymin,ymax],[xmin,xmax]]
		range1=(xmin,xmax,ymin,ymax)
	hh=scipy.histogram2d(y,x,range=range, bins=bins)
	plt.gca().set_xlabel(xtitle)
	plt.gca().set_ylabel(ytitle)

	plt.imshow(-hh[0],extent=range1, aspect='auto', interpolation='nearest')



def contour (z, x=None, y=None, xrange=None, yrange=None, zrange=None, xr=None, yr=None, zr=None, 
		title="", xtitle="", ytitle="", position=None, xlog=False, ylog=False, zlog=False, 
		xticklabel = "%1.1f", yticklabel = "%1.1f", zticklabel = "%1.1f", levels=None, 
		nlevels=256, c_color="black", cm ="jet", c_line="solid", c_levels=None, c_charsize=12.0, 
		c_thick=1, thick=1, font="monospace", weight="normal", charsize=14.0, bar=True, fill=True, 
		overplot=False, noerase=False):

# Adjust some built-in parameters to make plots nicer:		
#	matplotlib.rcParams['figure.facecolor']='white'
#	matplotlib.rcParams['text.usetex']=True
	matplotlib.rcParams['font.family']=font
	matplotlib.rcParams['font.size']=charsize
	matplotlib.rcParams['font.weight']=weight
	matplotlib.rcParams['axes.labelsize']=charsize	
	matplotlib.rcParams['xtick.labelsize']=charsize	
	matplotlib.rcParams['ytick.labelsize']=charsize			
	matplotlib.rcParams['lines.markeredgewidth']=thick-0.5	
	matplotlib.rcParams['lines.linewidth']=thick
	matplotlib.rcParams['lines.antialiased']=True	
	matplotlib.rcParams['xtick.minor.size']=7
	matplotlib.rcParams['ytick.minor.size']=7		
	matplotlib.rcParams['xtick.major.size']=12
	matplotlib.rcParams['ytick.major.size']=12		
#	matplotlib.rcParams['ps.papersize']=papersize


# Initialize x and y if these are not provided:
	if x==None or y==None:
		x, y = meshgrid(arange(z.shape[0]), arange(z.shape[1]))
		
# Define position of this plot:
	if not noerase:
		plt.gcf().clf()	
		if position!=None:
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
		z = log10(z)

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
	if levels==None:
		zmin = zrange[0]
		zmax = zrange[1]
		dz = zmax-zmin
		levels = zmin + dz/nlevels*arange(nlevels+1)	

# Setup frame thickness:
	plt.gca().frame.set_linewidth(thick) 

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
	xFormatter = ticker.FormatStrFormatter(xticklabel)
	yFormatter = ticker.FormatStrFormatter(yticklabel)
	plt.gca().xaxis.set_major_formatter(xFormatter)
	plt.gca().yaxis.set_major_formatter(yFormatter)

# Setup minor tickmarks:		
	xminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])
	yminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])	

	plt.gca().xaxis.set_minor_locator(xminorLocator)
	plt.gca().yaxis.set_minor_locator(yminorLocator)		
		
# Make a filled contour plot:		
	if fill:
		cset1=plt.gca().contourf(x, y, z, levels, cmap=plt.cm.get_cmap(cm, nlevels))

# Add contour lines:
	if c_levels == None:
		c_levels = levels[0:len(levels):int(nlevels/12)]
	cset2 = plt.gca().contour(x, y, z, c_levels, colors = c_color,linewidths=c_thick,hold='on')

# Do not display dashed contours for negative values:
	for c in cset2.collections:
   		c.set_linestyle(c_line)

# Add value labels on contour lines:
	zFormatter = ticker.FormatStrFormatter(zticklabel)
	cset3 = plt.gca().clabel(cset2, c_levels, inline=1, fmt=zticklabel, fontsize=c_charsize)
                
# Do we need a color bar?:                
	if bar:
		matplotlib.rcParams['ytick.labelsize']=c_charsize				
		plt.colorbar(cset1, ticks=[min(levels), max(levels)], shrink = 0.87, aspect=15, 
			format=zFormatter)
		
	if plt.isinteractive():
		plt.draw()
		
# Restore default matplotlib parameters:
	matplotlib.rcdefaults()


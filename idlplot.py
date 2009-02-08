import matplotlib.pyplot as plt
import numpy
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator

plt.ion()

#idlplot.plot(x,y,xrange=[0,39],yrange=[-1,10],ps=4,xtitle="X",\
#    color='black',position=[0.1,0.1,0.9,0.9], xlog=True)
#idlplot.oplot(x,2+y/10.,ps=3,xtitle="SX",color='blue')

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
	if overplot:
		tmp=oplot 
	else:
		tmp=plot
	xrange=[min,max]
	tmp(loc1,hh1,ps=0,color=color,xrange=xrange,yrange=yrange,
		xlog=xlog,ylog=ylog, **kw)
    
def plot (x,y, xrange=None, yrange=None, ps=0, thick=1, xtitle="", ytitle="",
		color='black', noerase=False, overplot=False,position=None, ylog=False,
		xlog=False):
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
		                    
	if ps==3:
		sym='.'
	elif ps==4: 
		sym='D'
	elif ps==7: 
		sym='+'
	elif ps==2: 
		sym='*'
	elif ps==6: 
		sym='s'
	else:
		sym='-'
	
	if xrange==None:
		xrange=[min(x),max(x)]
	if yrange==None:
		yrange=[min(y),max(y)]
#	xmajorLocator = MultipleLocator(1)
#	xmajorFormatter = FormatStrFormatter('%d')
#	ymajorLocator = MultipleLocator(1)
#	ymajorFormatter = FormatStrFormatter('%d')
	xminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])
	#MultipleLocator(.1)
#	yminorLocator = MultipleLocator(.5)
	yminorLocator = MaxNLocator(nbins=90, steps=[1, 2, 5, 10])

	plt.gca().xaxis.set_minor_locator(xminorLocator)
	plt.gca().yaxis.set_minor_locator(yminorLocator)		
		
	plt.gca().set_xlabel(xtitle)
	plt.gca().set_ylabel(ytitle)
	plt.gca().set_autoscale_on(False)
	if not overplot:
		plt.gca().axis(xrange+yrange)#,xautcoscale_on=False)
	plt.gca().plot(x,y,sym,linewidth=thick,color=color)
	if plt.isinteractive():
		plt.draw()
	
def oplot (x,y, **kw):

	plot (x,y, noerase=True, overplot=True, **kw)

#def ploterror (x,y, err **kw):
#	plot (x,y, **kw)
#	plt.gca().errorbar(x,y,err)
			

	
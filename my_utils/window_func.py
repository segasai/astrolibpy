import scipy,numpy

def window_func(x,y,func,xmin=None,xmax=None,nbin=100):
	"""This function does compute a user-supplied function 
	on the subsets of y grouped by values of x
	E.g. imagine that you have x from 0 to 100 and you want
	to know the mean values of y for 0<x<10, 10<x<2- etc..
	In that case you need to do
	window_func(x,y,lambda x: x.mean(),0,100,nbin=10)
	"""
	
	if xmin is None:
		xmin = x.min()
	if xmax is None:
		xmax = x.max()
	hh,loc=scipy.histogram(x,range=(xmin,xmax),bins=nbin)
	inds= numpy.digitize(x,loc)
	hh=hh*0.
	for i in range(nbin):
		hh[i]=func(y[inds==i])
	return (loc[1:]+loc[:-1])*0.5,hh
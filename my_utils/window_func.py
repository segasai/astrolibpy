import scipy,numpy

def window_func(x,y,func,xmin=None,xmax=None,nbin=100, empty=False ):
	"""This function does compute a user-supplied function 
	on the subsets of y grouped by values of x
	E.g. imagine that you have x from 0 to 100 and you want
	to know the mean values of y for 0<x<10, 10<x<2- etc..
	In that case you need to do
	xbin,funcy,nperbin=window_func(x,y,lambda x: x.mean(),0,100,nbin=10)
	where xbin is the array with the centers of bins, 
	funcy is the func -- evaluated for y where x is within the appropriate bin
	and nperbin is the number of points in the appropriate bin
	empty keyword is needed if you want to retrieve the function evaluation 
	in empty bins too, otherwise they aren't returned at all
	"""
	
	if xmin is None:
		xmin = x.min()
	if xmax is None:
		xmax = x.max()
	hh,loc=scipy.histogram(x,range=(xmin,xmax),bins=nbin)
	inds= numpy.digitize(x,loc)
	mask=numpy.zeros(len(hh),bool)
	retv=hh*0.	
	for i in range(nbin):
		cury=y[inds==i]
		if len(cury)>0 or empty:
			retv[i]=func(cury)
		mask[i]=len(cury)>0
	retx = (loc[1:]+loc[:-1])*0.5 # middle of the bin
	if empty:
		mask |= True
	return retx[mask], retv[mask],hh[mask]

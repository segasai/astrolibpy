import scipy,numpy,numpy as np

def window_func(x, y, func, xmin=None, xmax=None, nbin=100, empty=False,
			xlog=False):
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
	if xlog:
		xmin,xmax,x=[numpy.log10(tmp) for tmp in [xmin,xmax,x]]
	if (len(x)!=len(y) and x.ndim==1) or (x.shape!=y.shape and x.ndim>1):
		raise ValueError('Input arrays must have the same size')
	#hh,loc=scipy.histogram(x,range=(xmin,xmax),bins=nbin)
	inds = ((x-xmin)/float(xmax-xmin)*nbin).astype(int)
	mask = numpy.zeros(nbin, bool)
	retv = numpy.zeros(nbin)
	hh = numpy.zeros(nbin,int)

	for i in range(nbin):
		cury=y[inds==i]
		hh[i]=len(cury)
		mask[i]=len(cury)>0
		if len(cury)>0 or empty:
			retv[i]=func(cury)
	retx = xmin+(xmax-xmin)*1./nbin*(0.5+numpy.arange(nbin))
	if xlog:
		retx = 10**retx
	if empty:
		mask |= True
	return retx[mask], retv[mask],hh[mask]


def window_func2d(x,y,z,func,xmin=None,xmax=None,ymin=None,ymax=None,nbins=[10,10]):
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
	if ymin is None:
		ymin = y.min()
	if ymax is None:
		ymax = y.max()

	hh, locx, locy = scipy.histogram2d(x, y,
			range=((xmin,xmax), (ymin,ymax)), bins=nbins)
	xinds = numpy.digitize(x, locx) - 1
	yinds = numpy.digitize(y, locy) - 1
	mask = numpy.zeros(hh.shape, bool)
	retv = hh  * 0.	
	subind = (xinds >= 0) & (yinds >= 0) & (xinds < nbins[0]) & (yinds < nbins[1])
	xinds, yinds, z1 = [_[subind] for _ in [xinds,yinds,z]]
	valind = yinds * (nbins[0]) + xinds
	sortind = np.argsort(valind)
	valind = valind[sortind]
	poss = np.where(np.diff(valind) > 0)[0]
	z1 = z1[sortind]
	for  i in range(len(poss) + 1):
		if i == 0:
			left = 0
			right = poss[0] + 1
		elif i == len(poss):
			left = poss[-1] + 1
			right = len(valind)
		else:
			left = poss[i-1] + 1
			right = poss[i] + 1
		curval = valind[left]
		retv[curval % (nbins[0]), curval / (nbins[0])] = func(z1[left:right])
	return retv, locx[0], locx[-1], locy[0], locy[-1], hh
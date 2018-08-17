import idlplot
import numpy as np

""" 
This module is a set wrappers around idlplot designed to make 
the plots of subsets of the data: e.g.
plot(x,y,ind=ind)
instead of plot(x[ind],y[ind])
  
"""
def tvhist2d(a,b,*args,**kw):
	ind = kw.get('ind')
	
	if ind is None:
		return idlplot.tvhist2d(a,b,*args,**kw)
	else:
		weights = kw.get('weights')
		if weights is not None:
			kw['weights']=kw['weights'][ind]
		del kw['ind']
		return idlplot.tvhist2d(a[ind],b[ind],*args,**kw)

def plothist(a,*args,**kw):
	ind = kw.get('ind')
	
	if ind is None:
		ret=idlplot.plothist(a,*args,**kw)
	else:
		weights = kw.get('weights')
		if weights is not None:
			if not np.isscalar(kw['weights']):
				kw['weights']=kw['weights'][ind]
		del kw['ind']
		ret=idlplot.plothist(a[ind],*args,**kw)
	return ret

def plot(a, b=None, **kw):
	ind = kw.get('ind')
	
	if ind is None:
		idlplot.plot(a,b,**kw)
	else:
		del kw['ind']
		if b is not None:
			idlplot.plot(a[ind], b[ind], **kw)
		else:
			idlplot.plot(a[ind], None, **kw)
			

def oplot(a, b=None, **kw):
	ind = kw.get('ind')
	
	if ind is None:
		idlplot.oplot(a, b, **kw)
	else:
		del kw['ind']
		if b is not None:
			idlplot.oplot(a[ind], b[ind], **kw)
		else:
			idlplot.oplot(a[ind], **kw)

def errorfixer(var, ind):
	var = np.asarray(var)
	if var.ndim==2 and var.shape[0]==2:
		var1 = [var[0][ind],var[1][ind]]
	else:
		var1 = var[ind]
	return var1

def ploterror(a,b,c,*args,**kw):
	ind = kw.get('ind')
	
	if ind is None:
		idlplot.ploterror(a,b,c,*args,**kw)
	else:
		del kw['ind']
		l = len(args)
		args1=[None]*l
		c1 = errorfixer(c, ind)

		for i in range(l):
			args1[i] = errorfixer(args[i],ind)
		idlplot.ploterror(a[ind],b[ind],c1,*args1,**kw)

tvhist2d.__doc__ = idlplot.tvhist2d.__doc__
plot.__doc__ = idlplot.plot.__doc__
oplot.__doc__ = idlplot.oplot.__doc__
ploterror.__doc__ = idlplot.ploterror.__doc__
plothist.__doc__ = idlplot.plothist.__doc__

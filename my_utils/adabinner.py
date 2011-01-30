import numpy.random, numpy,quick_hist, scipy.stats, random
from idlsave import idlsave

def doit(x, y, hi1=None, hi2=None, thresh=None):
	hhs ={}
	for curhi in range(hi1,hi2+1):
		hh = quick_hist.quick_hist((x, y), range=((0, 1), (0, 1)),
								nbins=[2**curhi,2**curhi])
		hhs[curhi]=hh

	hh0 = hhs[hi2] * 1 # accumulator of the result
	area = hh0 * 0

	two = 2
	poiss = scipy.stats.poisson(thresh)

	DeepenOrNot = lambda x: random.random() < poiss.cdf(x)
	#DeepenOrNot = lambda x: x>thresh

	def doitit(it, i, j):
		curhhs = hhs[it]
		if it==hi2:
			hh[i:i+2, j:j+2] = curhhs[i:i+2, j:j+2]
			area[i:i+2, j:j+2] = 2**(-it)
		else:
			for ii in range(i, i+2):
				for jj in range(j, j+2):
					curval = curhhs[ii,jj]
					if DeepenOrNot(curval):
						doitit(it+1, ii*two, jj*two)
					else:
						dx=2**(hi2-it)
						hh0[ii*dx:(ii+1)*dx, jj*dx:(jj+1)*dx]=curval
						area[ii*dx:(ii+1)*dx, jj*dx:(jj+1)*dx] = 2**(-it)
	n1 = 2**hi1
	dn = 2**(hi2-hi1)

	for i in range(n1):
		for j in range(n1):
				if DeepenOrNot(hhs[hi1][i,j]):
					doitit(hi1+1,i*two,j*two)
				else:
					hh0[i*dn:(i+1)*dn, j*dn:(j+1)*dn] = hhs[hi1][i,j]
					area[i*dn:(i+1)*dn, j*dn:(j+1)*dn] = 2**(-hi1)					
	return hh0*1./(2**hi1*area)**2							

def hist2d(x, y, xmin=None, xmax=None, ymin=None, ymax=None, hi=[2,10],
			thresh=30):

	xmin = x.min() if xmin is None else xmin
	ymin = y.min() if ymin is None else ymin
	xmax = x.max() if xmax is None else xmax
	ymax = y.max() if ymax is None else ymax

	xmod = (x-xmin)/(xmax-xmin)
	ymod = (y-ymin)/(ymax-ymin)

	ind = (xmod>=0)&(xmod<=1)&(ymod>=0)&(ymod<=1)
	res = doit(xmod[ind], ymod[ind], hi1=hi[0], hi2=hi[1], thresh=thresh)	
	return res


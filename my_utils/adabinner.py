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


import numpy.random, numpy,quick_hist, scipy.stats, random


__doc__="""
Adaptive binner module
Functions:
hist:
	1D adaptive histogram
hist2d:
	2D adaptive histogram

"""
def __doit2d(x, y, hi1=None, hi2=None, thresh=None):
	hhs ={}
	for curhi in range(hi1,hi2+1):
		hh = quick_hist.quick_hist((x, y), range=((0, 1), (0, 1)),
								nbins=[2**curhi,2**curhi])
		hhs[curhi]=hh
	pixcen = []
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
			pixcen.append((i+.5,j+.5))
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
						pixcen.append((ii*dx+dx/2.+0.5,jj*dx+dx/2.+0.5))
	n1 = 2**hi1
	dn = 2**(hi2-hi1)

	for i in range(n1):
		for j in range(n1):
				if DeepenOrNot(hhs[hi1][i,j]):
					doitit(hi1+1,i*two,j*two)
				else:
					hh0[i*dn:(i+1)*dn, j*dn:(j+1)*dn] = hhs[hi1][i,j]
					area[i*dn:(i+1)*dn, j*dn:(j+1)*dn] = 2**(-hi1)					
					pixcen.append((i*dn+dn/2.+0.5,j*dn+dn/2.+0.5))
	area = (2**hi1*area)**2 # in smallest pixels squared
	return hh0*1./area,pixcen,area


def __doit1d(x, hi1=None, hi2=None, thresh=None):
	hhs ={}
	for curhi in range(hi1,hi2+1):
		hh = quick_hist.quick_hist((x,), range=((0, 1),),
								nbins=[2**curhi])
		hhs[curhi]=hh

	hh0 = hhs[hi2] * 1 # accumulator of the result
	area = hh0 * 0

	two = 2
	poiss = scipy.stats.poisson(thresh)

	DeepenOrNot = lambda x: random.random() < poiss.cdf(x)
	#DeepenOrNot = lambda x: x>thresh

	def doitit(it, i):
		curhhs = hhs[it]
		if it==hi2:
			hh[i:i+2] = curhhs[i:i+2]
			area[i:i+2] = 2**(-it)
		else:
			for ii in range(i, i+2):
				curval = curhhs[ii]
				if DeepenOrNot(curval):
					doitit(it+1, ii*two)
				else:
					dx=2**(hi2-it)
					hh0[ii*dx:(ii+1)*dx]=curval
					area[ii*dx:(ii+1)*dx] = 2**(-it)
	n1 = 2**hi1
	dn = 2**(hi2-hi1)

	for i in range(n1):
			if DeepenOrNot(hhs[hi1][i]):
				doitit(hi1+1,i*two)
			else:
				hh0[i*dn:(i+1)*dn] = hhs[hi1][i]
				area[i*dn:(i+1)*dn] = 2**(-hi1)
	return hh0*1./(2**hi1*area)							


def hist2d(x, y, xmin=None, xmax=None, ymin=None, ymax=None, hi=[2,10],
			thresh=30, full_output=False):
	"""
	This function does the 2D histogram with adaptive binning 
	Example:
	>> hist2d(xs,ys, hi=[3,6], thresh=30)
	>> hh,xloc,yloc,pixcen,area = hist2d(xra,dec,full_output=True,thresh=100)                
	Keyword parameters:
	------------------
	hi
		the list of two integer values: they describe how coarse the
		largest bin and how fine is the smallest bin, e.g. 
		[2,5] means the largest possible bin has a size 
		of 1/2**2 of the your dataset and the smallest bin has a
		size of 1/2**5
	thresh
		the minimum number of points within a bin allowed (
		if the number is smaller than the threshold then further
		decreasing of the bin size is not allowed by the algorithm)
	xmin,xmax,ymin,ymax
		x-ranges and y-ranges. If not specified, they are determined
		from the x.min(),x.max(),y.min(),y.max()
	full_output
		Boolean controlling whether to output just just the histogram(full_output=False)
		or output the tuple with the histogram, grid-centers in x, grid-centers in y,
		pixel-centers, and area
	"""
	xmin = x.min() if xmin is None else xmin
	ymin = y.min() if ymin is None else ymin
	xmax = x.max() if xmax is None else xmax
	ymax = y.max() if ymax is None else ymax

	xmod = (x-xmin)/(xmax-xmin)
	ymod = (y-ymin)/(ymax-ymin)

	ind = (xmod>=0)&(xmod<=1)&(ymod>=0)&(ymod<=1)
	hh,pixcen,area = __doit2d(xmod[ind], ymod[ind], hi1=hi[0], hi2=hi[1], thresh=thresh)	
	xloc = numpy.linspace(xmin,xmax,hh.shape[0],False)
	yloc = numpy.linspace(ymin,ymax,hh.shape[0],False)
	xloc += 0.5*(xloc[1]-xloc[0])
	yloc += 0.5*(yloc[1]-yloc[0])
	
	if full_output:
		out = hh,xloc,yloc,pixcen,area
	else:
		out = hh
	return out

def hist(x, xmin=None, xmax=None, hi=[2,10], thresh=30):
	"""
	This function does the 1D histogram with adaptive binning 
	Example:
	>> loc, hh = hist(xs, hi=[3,6], thresh=30)
	
	Keyword parameters:
	------------------
	hi
		the list of two integer values: they describe how coarse the
		largest bin and how fine is the smallest bin, e.g. 
		[2,5] means the largest possible bin has a size 
		of 1/2**2 of the your dataset and the smallest bin has a
		size of 1/2**5
	thresh
		the minimum number of points within a bin allowed (
		if the number is smaller than the threshold then further
		decreasing of the bin size is not allowed by the algorithm)
	xmin,xmax
		x-range. If not specified, they are determined
		from the x.min(),x.max()
	--------------------------------------------	
	Returns:
		the histogram and the bin edges vector
	"""
	xmin = x.min() if xmin is None else xmin
	xmax = x.max() if xmax is None else xmax

	xmod = (x-xmin)/(xmax-xmin)

	ind = (xmod>=0)&(xmod<=1)
	hh = __doit1d(xmod[ind], hi1=hi[0], hi2=hi[1], thresh=thresh)	
	loc = numpy.linspace(xmin,xmax,len(hh)+1,True)
	return hh, loc


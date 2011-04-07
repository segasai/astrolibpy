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


import numpy
import scipy.weave, scipy.weave.converters

def quick_hist(arrs, range=None, nbins=None, weights=None):
	"""
	arr -- tuple of N-arrays
	range -- list of tuples of ranges
	nbins -- list of lengths	
	"""
	from __builtin__ import range as xrange
	nd = len(arrs)
	if len(nbins)!=nd:
		raise ValueError('The array of nbins MUST have the same length as the number of input data vectors')
	if len(range)!=nd:
		raise ValueError('The array of ranges MUST have the same length as the number of input data vectors')

	nx = len(arrs[0])
	for curarr in arrs:
		if len(curarr)!=nx:
			raise ValueError('All the input arrays MUST have the same length!')
	if weights is not None:
		if len(weights)!=nx:
			raise ValueError('The weights array MUST have the same length as the input arrays')
	
	poss = numpy.zeros((nx,), dtype=numpy.int64)

	ind = numpy.ones_like(arrs[0]).astype(bool)
	nbins_rev = nbins + []
	nbins_rev.reverse()
	mults = (reduce(lambda x, y: x + [y * x[-1]], nbins_rev, [1]))[:-1]
 	mults.reverse()	
	for i in xrange(nd):
		cur_arr = arrs[i]
		cur_range0 = range[i][0]
		cur_range1 = range[i][1]
		cur_nbins = nbins[i]
		cur_mult = mults[i]
		code1 = """
		int i, cur_pos;
		double curfac = cur_nbins * 1./ (cur_range1-cur_range0);
		for (i=0; i<nx; i++)
		{
			cur_pos = floor( ( cur_arr(i)-cur_range0) * curfac);

			if ((cur_pos>=0 ) && (cur_pos<cur_nbins))
			{
				poss(i)=poss(i)+cur_pos*cur_mult;
			}
			else
			{
				ind(i)=false;
			}
		}"""
		try:
			if weights is None:
				scipy.weave.inline(code1, ['cur_range0', 'cur_range1', 'cur_nbins','poss','ind','nx','cur_arr','cur_mult'],
					type_converters=scipy.weave.converters.blitz)
		except:
			cur_pos = (cur_arr - cur_range0) * (cur_nbins * 
									1. / (cur_range1 - cur_range0))
			cur_pos = numpy.floor(cur_pos).astype(numpy.int64)
			ind &= ((cur_pos >= 0) & ( cur_pos < cur_nbins))
			poss += cur_pos * cur_mult

	poss = poss[ind]
	newlen = len(poss)
	if weights is None:
		weights_str = '1'
	else:
		weightsind = weights[ind]
		weights_str = 'weightsind(i)'
	del ind
	res = numpy.zeros(numpy.array(nbins, dtype=numpy.int64).prod())

	code = """
	int i;
	for (i=0; i<newlen; i++)
	{
		res(poss(i)) = res(poss(i)) + %s;
	}"""%weights_str
	try:
		if weights is None:
			scipy.weave.inline(code, ['res', 'poss', 'newlen'],
				type_converters=scipy.weave.converters.blitz)
		else:
			scipy.weave.inline(code, ['res', 'poss', 'newlen','weightsind'],
				type_converters=scipy.weave.converters.blitz)			
	except Exception:
		print "Sorry the compiled version didn't work :("
		if weights is None:
			for i in xrange(len(poss)):
				res[poss[i]]+=1
		else:
			for i in xrange(len(poss)):
				res[poss[i]]+=weights[i]
	
	return res.reshape(nbins)
		
		
# Copyright (C) 2009-2010 Sergey Koposov
# This file is part of astrolibpy
#
#	 astrolibpy is free software: you can redistribute it and/or modify
#	 it under the terms of the GNU General Public License as published by
#	 the Free Software Foundation, either version 3 of the License, or
#	 (at your option) any later version.
#
#	astrolibpy is distributed in the hope that it will be useful,
#	 but WITHOUT ANY WARRANTY; without even the implied warranty of
#	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	 GNU General Public License for more details.
#
#	 You should have received a copy of the GNU General Public License
#	 along with astrolibpy.	 If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function
from functools import reduce
import warnings
import numpy as np
import hashlib
import importlib
import os

def _getcode():
	cdef="""
	void hister(
		double cur_range0,
		double cur_range1,
		int cur_nbins,
		int cur_mult,
		double *cur_arr,
		int64_t *poss,
		int64_t *ind,
		int nx);

	void adder_now(
		int64_t* res,
		int64_t* poss,
		int newlen);

	void adder_wei(
		double* res,
		int64_t* poss,
		double* weights,
		int newlen);

	"""								

	src="""
	static void hister(
		double cur_range0,
		double cur_range1,
		int cur_nbins,
		int cur_mult,
		double *cur_arr,
		int64_t *poss,
		int64_t *ind,
		int nx)
		{
			int64_t i, cur_pos;
			double curfac = cur_nbins * 1./ (cur_range1- cur_range0);
			for(i=0;i<nx;i++)
			{
				cur_pos = (int)floor( ( cur_arr[i]-cur_range0) * curfac);
				if ((cur_pos>=0 ) && (cur_pos<cur_nbins))
				{
					poss[i] = poss[i] + cur_pos * cur_mult;
				}
				else
				{
					ind[i]=0;
				}    		 
				 
			}
		}

	static void adder_now(
		int64_t* res,
		int64_t* poss,
		int newlen)
	{
		int i;
		for (i=0; i<newlen ;i++)
		{
			res[poss[i]] = res[poss[i]] + 1;
		}
	}

	static void adder_wei(
		double* res,
		int64_t* poss,
		double* weights,
		int newlen)
	{
		int i;
		for (i=0; i<newlen ;i++)
		{
			res[poss[i]] = res[poss[i]]+weights[i];
		}
	}
		
	"""
	return cdef, src

def _hasher(cdef, src):
	md5 = hashlib.md5((cdef+src).encode('utf-8')).hexdigest()[:6]
	return md5

def _getlibname():
	hash = _hasher(*_getcode())
	return '_quick_hist_'+hash

def _buildlib():
	cdef, src = _getcode()
	libname = _getlibname()
	from cffi import FFI
	ffibuilder = FFI()
	ffibuilder.cdef(cdef)
	ffibuilder.set_source(libname, src)
	dir = os.path.abspath(os.path.dirname(__file__))
	ffibuilder.compile(tmpdir=dir,verbose=True)


def quick_hist(arrs, range=None, nbins=None, weights=None, getPos=False):
	"""
	N-dimensional histogram routine.
	Example:
	> xs=np.random.uniform(size=100); ys= np.random.uniform(size=100)
	> hh = quick_hist((xs,ys), range=[(0,1),(0,1)], nbins=[20,10])
	Arguments:
			arr -- tuple of N-arrays
			range -- list of tuples of ranges
			nbins -- list of numbers of bins 
	Keywords:
			weights -- weighting for the histogram
			getPos -- return the 1D vector of the positions within the histogram
									(-1 if the point is outside the range)
	"""
	nd = len(arrs)
	if range is None:
		range = []
		for i in np.arange(nd):
			range.append((arrs[0].min(), arrs[0].max()))
	if nbins is None:
		nbins = [10] * nd

	if len(nbins) != nd:
		raise ValueError(
			'The array of nbins MUST have the same length as the number of input data vectors')
	if len(range) != nd:
		raise ValueError(
			'The array of ranges MUST have the same length as the number of input data vectors')

	nx = arrs[0].size
	for curarr in arrs:
		if (curarr.size) != nx:
			raise ValueError('All the input arrays MUST have the same length!')
	if weights is not None:
		if (weights.size) != nx:
			raise ValueError(
				'The weights array MUST have the same length as the input arrays')
	# convert all the bins into integers
	nbins = [int(_tmp) for _tmp in nbins]

	poss = np.zeros((nx,), dtype=np.int64)
	ind = np.ones((nx,), dtype=np.int64)
	nbins_rev = nbins + []
	nbins_rev.reverse()
	mults = (reduce(lambda x, y: x + [y * x[-1]], nbins_rev, [1]))[:-1]
	mults.reverse()
	assert(poss.flags['C_CONTIGUOUS'])
	assert(ind.flags['C_CONTIGUOUS'])
	slow = False
	modname = _getlibname()

	for i in np.arange(nd):
		cur_arr = np.ascontiguousarray(arrs[i], dtype=np.float64)
		cur_range0 = float(range[i][0])
		cur_range1 = float(range[i][1])
		cur_nbins = int(nbins[i])
		cur_mult = int(mults[i])
		try:
			try:
				M = importlib.import_module(modname)
				#from _quick_hist import ffi, lib
			except:
				_buildlib()
				M = importlib.import_module(modname)
			poss_ffi = M.ffi.cast('int64_t *', M.ffi.from_buffer(poss))
			ind_ffi = M.ffi.cast('int64_t *', M.ffi.from_buffer(ind))
			cur_arr_ffi = M.ffi.cast('double *', M.ffi.from_buffer(cur_arr))

			M.lib.hister(cur_range0, cur_range1, cur_nbins,
					   cur_mult, cur_arr_ffi,
					   poss_ffi, ind_ffi, nx)
		except:
			warnings.warn("""Sorry the compiled version didn't work :(, executing a slower Python-only version.""")
			slow = True
			cur_pos = (cur_arr - cur_range0) * (cur_nbins *
												1. / (cur_range1 - cur_range0))
			cur_pos = np.floor(cur_pos).astype(np.int64)
			ind = ind * ((cur_pos >= 0) & (cur_pos < cur_nbins))
			poss += cur_pos * cur_mult

	ind = ind.astype(bool)

	poss = poss[ind]
	newlen = len(poss)


	nret = np.array(nbins, dtype=np.int64).prod()

	if weights is not None:
		res = np.zeros(nret, dtype=np.float64)
		weightsind = np.ascontiguousarray(weights[ind], dtype=np.float64)
	else:
		res = np.zeros(nret, dtype=np.int64)

	if not getPos:
		del ind

	assert(res.flags['C_CONTIGUOUS'])
	assert(poss.flags['C_CONTIGUOUS'])

	if not slow:
		poss_ffi = M.ffi.cast('int64_t *', M.ffi.from_buffer(poss))

		if weights is not None:
			weightsind_ffi = M.ffi.cast('double *', M.ffi.from_buffer(weightsind))
			res_ffi = M.ffi.cast('double *', M.ffi.from_buffer(res))
			M.lib.adder_wei(res_ffi, poss_ffi, weightsind_ffi, newlen)
		else:
			res_ffi = M.ffi.cast('int64_t *', M.ffi.from_buffer(res))
			M.lib.adder_now(res_ffi, poss_ffi, newlen)
	else:
		if weights is None:
			for i in np.arange(len(poss)):
				res[poss[i]] += 1
		else:
			for i in np.arange(len(poss)):
				res[poss[i]] += weightsind[i]
	if not getPos:
		return res.reshape(nbins)
	else:
		H = np.zeros(len(ind), dtype=np.int64) - 1
		H[ind] = poss
		return res.reshape(nbins), H

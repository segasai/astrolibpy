# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 14:00:40 2013

@author: koposov
"""
import numpy as np
import workerpool 

class derivator:
	# class which compute finite difference derivatives using multiple cpus
	def __init__(self, func, eps1=1e-4, eps2=1e-8, nthreads=1 ,kw={}):
		self.func=func
		if nthreads > 1:
			self.pool = workerpool.pool(func, nthreads=nthreads, kw=kw)
		self.eps1 = eps1
		self.eps2 = eps2

	def __call__(self,p0):
		p = np.asarray(p0)
		eps = np.maximum(np.abs(p * self.eps1), self.eps2)
		if self.pool is not None:
			self.pool.apply_async(0, p)
			for i in range(len(p)):
				curp = p.copy()
				curp[i] += eps[i]			
				self.pool.apply_async(i + 1, curp)
			res0 = self.pool.get(0)
			res = np.zeros(len(p))
			for i in range(len(p)):
				res[i] = self.pool.get(i+1)
		else:
			res0 = self.func(p)
			for i in range(len(p)):
				curp = p.copy()
				curp[i] += eps[i]
				res[i] = self.func(curp)
		ret = (res - res0) / eps
		return ret

	def __del__(self):
		if self.pool is not None:
			del self.pool
		del self.func
			

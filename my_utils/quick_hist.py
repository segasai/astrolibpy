import numpy
import scipy.weave, scipy.weave.converters

def quick_hist(arrs, range=None, nbins=None):
	from __builtin__ import range as xrange
	nd = len(arrs)
	nx = len(arrs[0])
	poss = numpy.zeros((nx,), dtype=numpy.int64)
	ind = numpy.ones_like(arrs[0]).astype(bool)
	nbins_rev = nbins + []
	nbins_rev.reverse()
	mults = (reduce(lambda x, y: x + [y * x[-1]], nbins_rev, [1]))[:-1]
 	mults.reverse()	
	for i in xrange(nd):
		cur_arr = arrs[i]
		cur_range = range[i]
		cur_nbins = nbins[i]
		cur_pos = (cur_arr - cur_range[0]) * (cur_nbins * 
								1. / (cur_range[1] - cur_range[0]))
		cur_pos = cur_pos.astype(numpy.int64)
		ind = (ind) & (cur_pos >= 0) & ( cur_pos < cur_nbins)
		poss = poss + cur_pos * mults[i]
	poss = poss[ind]
	res = numpy.zeros(numpy.array(nbins, dtype=numpy.int64).prod())
	newlen = len(poss)
	code = """
	int i;
	for (i=0; i<newlen; i++)
	{
		res(poss(i)) = res(poss(i)) + 1;
	}"""
	try:
		scipy.weave.inline(code, ['res', 'poss', 'newlen'],
			type_converters=scipy.weave.converters.blitz)
	except Exception:
		print "Sorry the compiled version didn't work :("
		for i in xrange(len(poss)):
			res[poss[i]]+=1
	
	return res.reshape(nbins)
		
		
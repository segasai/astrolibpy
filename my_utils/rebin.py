import scipy.signal,numpy as np

def rebin(arr, shape1):
	# rebins the array too the new shape by summing the elements
	# it is assumed that in each dimension the resizing is by an integer
	# factor
	shape0 = arr.shape 
	rats = np.array([_//__ for _, __ in zip (shape0,shape1)])
	for r, s0, s1 in zip(rats,shape0,shape1):
		if r*s1!=s0:
			raise Exception('Wrong array sizes')
	conv=scipy.signal.convolve(arr,
		np.ones(dtype=np.int64,shape=rats),'valid')
	ret = conv[tuple([slice(0,None,_) for _ in rats])]
	return ret
	
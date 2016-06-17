import multiprocessing as mp, numpy as np
import signal

def doit(func, queuein, queout, args, kw):	
	signal.signal(signal.SIGINT, signal.SIG_IGN)
	while True:
		inp = queuein.get()
		if inp is None:
			break
		i, p = inp
		res = func(p, *args, **kw)
		queout.put((i, res))

class pool:
	def __init__(self, func, args=None, kw={}, nthreads=1, mpmap=False):
		# create a processing pool that has the function that needs to 
		# be called stored
		# that reduces the overhead of constantly pickling the function
		# mpmap=True flag enables the multiprocessing.Pool.map syntax of map
		# e.g. map(func, list) 
		# the first argument is then ignored
		self.nthreads = nthreads
		self.queuein = mp.Queue()
		self.queueout = mp.Queue()
		self.mpmap = mpmap
		self.procs = [mp.Process(target=doit, name='xpool_%d'%i,
				args=(func, self.queuein, self.queueout, args, kw))
					for i in range(nthreads) ]
		[_.start() for _ in self.procs]
		self.stage = {}

	def apply_async(self, i, p):
		self.queuein.put((i, p))

	def get(self, i):
		if i in self.stage:
			ret = self.stage[i]
			del self.stage[i]
			return ret
		else:
			while True:
				iret, ret = self.queueout.get()
				if iret != i:
					self.stage[iret] = ret
				else:
					return ret

	def get_any(self):
		if len(self.stage) != 0:
			i, ret = self.stage.popitem()
			return i, ret
		else:
			iret, ret = self.queueout.get()
			return iret, ret

	def map(self, *args):
		if self.mpmap:
			params=args[1]
		else:
			params=args[0]
		for i, p in enumerate(params):
			self.queuein.put((i, p))
		ret = [None] * len(params)
		for i in range(len(params)):
			resi, res = self.queueout.get()
			ret[resi] = res
		return ret

	def __del__(self):
		self.join()
			
	def join(self):
		for p in self.procs:
			self.queuein.put(None)
		for p in self.procs:
			p.join()
		del self.queuein
		del self.queueout

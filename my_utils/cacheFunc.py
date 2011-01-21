from collections import deque


class cacheFunc:
	""" class to cache function calls
	E.g.
	>> def func(a,b=3):
	>>	# do something 
	>>	print 'Called', a,b
	>>	return a+b
		
	>> funcC = cacheFunc(func)
	>> print funcC(2,b=3)
	Called 2 3
	5
	>> print funcC(3,b=4)
	Called 3 4
	7
	>> print funcC(2,b=3)
	5
	It also can be used as decorator:
	@cacheFunc
	def myfunc(a,b,c,d=3,e=4):
	....
	"""
	def __init__(self, func,maxn=10):
		self.hash={}
		self.l=deque()
		self.maxn=maxn
		self.func=func
	
	def __call__(self,*args, **kw):
		key =  args+tuple(kw.iteritems())
		if self.hash.get(key) is not None:
			pass
		else:
			self.hash[key]=self.func(*args, **kw)
			self.l.append(key)
			if len(self.l)>self.maxn:
				delarg = self.l.popleft()
				del self.hash[delarg]
				del delarg
		return self.hash[key]
								
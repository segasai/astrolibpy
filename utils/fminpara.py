from scipy.optimize import OptimizeResult
import numpy as np
import workerpool

_status_message = {'success': 'Optimization terminated successfully.',
                   'maxfev': 'Maximum number of function evaluations has '
                              'been exceeded.',
                   'maxiter': 'Maximum number of iterations has been '
                              'exceeded.',
                   'pr_loss': 'Desired error not necessarily achieved due '
                              'to precision loss.'}

def fmin_para(func, x0, xtol=1e-4, ftol=1e-4, maxiter=None, maxfun=None,
		 full_output=0, disp=1, retall=0, callback=None, nthreads=1, nverts=None,
			nactive=None):
	"""
	Minimize a function using the downhill simplex algorithm.
	This code does it in parallel fashion (modified by Sergey Koposov)
	
	This algorithm only uses function values, not derivatives or second
	derivatives.

	Parameters
	----------
	func : callable func(x,*args)
		The objective function to be minimized.
	x0 : ndarray
		Initial guess.
	args : tuple, optional
		Extra arguments passed to func, i.e. ``f(x,*args)``.
	callback : callable, optional
		Called after each iteration, as callback(xk), where xk is the
		current parameter vector.
	xtol : float, optional
		Relative error in xopt acceptable for convergence.
	ftol : number, optional
		Relative error in func(xopt) acceptable for convergence.
	maxiter : int, optional
		Maximum number of iterations to perform.
	maxfun : number, optional
		Maximum number of function evaluations to make.
	full_output : bool, optional
		Set to True if fopt and warnflag outputs are desired.
	disp : bool, optional
		Set to True to print convergence messages.
	retall : bool, optional
		Set to True to return list of solutions at each iteration.

	Returns
	-------
	xopt : ndarray
		Parameter that minimizes function.
	fopt : float
		Value of function at minimum: ``fopt = func(xopt)``.
	iter : int
		Number of iterations performed.
	funcalls : int
		Number of function calls made.
	warnflag : int
		1 : Maximum number of function evaluations made.
		2 : Maximum number of iterations reached.
	allvecs : list
		Solution at each iteration.

	See also
	--------
	minimize: Interface to minimization algorithms for multivariate
		functions. See the 'Nelder-Mead' `method` in particular.

	Notes
	-----
	Uses a Nelder-Mead simplex algorithm to find the minimum of function of
	one or more variables.

	This algorithm has a long history of successful use in applications.
	But it will usually be slower than an algorithm that uses first or
	second derivative information. In practice it can have poor
	performance in high-dimensional problems and is not robust to
	minimizing complicated functions. Additionally, there currently is no
	complete theory describing when the algorithm will successfully
	converge to the minimum, or how fast it will if it does.

	References
	----------
	.. [1] Nelder, J.A. and Mead, R. (1965), "A simplex method for function
		   minimization", The Computer Journal, 7, pp. 308-313

	.. [2] Wright, M.H. (1996), "Direct Search Methods: Once Scorned, Now
		   Respectable", in Numerical Analysis 1995, Proceedings of the
		   1995 Dundee Biennial Conference in Numerical Analysis, D.F.
		   Griffiths and G.A. Watson (Eds.), Addison Wesley Longman,
		   Harlow, UK, pp. 191-208.

	"""
	opts = {'xtol': xtol,
			'ftol': ftol,
			'maxiter': maxiter,
			'maxfev': maxfun,
			'disp': disp,
			'return_all': retall,
			'nthreads':nthreads,
			'nverts':nverts,
			'nactive':nactive}

	res = _minimize_neldermead_para(func, x0, callback=callback,  **opts)

	if full_output:
		retlist = res['x'], res['fun'], res['nit'], res['nfev'], res['status']
		if retall:
			retlist += (res['allvecs'], )
		return retlist
	else:
		if retall:
			return res['x'], res['allvecs']
		else:
			return res['x']

def _minimize_neldermead_para(func, x0, callback=None,
						 xtol=1e-4, ftol=1e-4, maxiter=None, maxfev=None,
						 disp=False, return_all=False,pool=None, nthreads=1,
						nverts=None, nactive=None):
	maxfun = maxfev
	retall = return_all
	allvecs = []

	x0 = np.asfarray(x0).flatten()
	N = len(x0)

	fcalls = [0]
	allcalls = []
	if nverts is None:
		nverts = (N+1)*2
	assert(nverts>=(N+1))

	if nthreads>=nverts:
		nthreads = nverts
	pool = workerpool.pool(func, nthreads)
	if nactive is None:
		nactive = nthreads # number of active points
	if nactive >= nverts:
		nactive = nverts-1
	def applicator(i, x):
		#print x
		pool.apply_async(i, x)
		fcalls[0] +=1
		allcalls.append(x.tolist())

	rank = len(x0.shape)
	if not -1 < rank < 2:
		raise ValueError("Initial guess must be a scalar or rank-1 sequence.")
	if maxiter is None:
		maxiter = N * 200
	if maxfun is None:
		maxfun = N * 200

	rho = 1
	chi = 2
	psi = 0.5
	sigma = 0.5
	nonzdelt = 0.05
	zdelt = 0.00025
	excess = nverts - N -1 
	one2np1 = list(range(1, nverts))

	if rank == 0:
		sim = np.zeros((nverts,), dtype=x0.dtype)
	else:
		sim = np.zeros((nverts, N), dtype=x0.dtype)
	fsim = np.zeros((nverts,), float)
	
	applicator(0, x0)
	sim[0] = x0
	if False:
		for k in range(0, N):
			y = np.array(x0, copy=True)
			if y[k] != 0:
				y[k] = (1 + nonzdelt) * y[k]
			else:
				y[k] = zdelt
			sim[k + 1] = y
			applicator(k + 1, y)
			
		nover = int(np.ceil(excess *1. / N))
		for k in range(N, nverts-1):
			pos1 = (k-N) / nover
			pos2 = (k-N) % nover
			#print k,pos1,pos2
			y = (sim[1+pos1]-x0) * (pos2+1)*1./(nover+1)+x0
			#y = np.array(x0, copy=True)	
			#fac= np.random.uniform(0,zdelt,size=N)
			#y = (x0*(1+fac))*(x0!=0).astype(int)+(fac)*(x0==0).astype(int)
			sim[k + 1] = y
			applicator(k + 1, y)
	state=np.random.get_state()
	np.random.seed(1)

	for k in range(0, nverts-1):
		y = np.array(x0, copy=True)
		rand=np.random.normal(0,1, len(x0))
		rand=rand/(rand**2).sum()**.5+1./len(x0)**.5
		y = (x0)*(1+rand*nonzdelt)*(x0!=0).astype(int)+(x0==0).astype(int)*(rand)*nonzdelt
		sim[k + 1] = y
		applicator(k + 1, y)
	np.random.set_state(state)
	for k in range(0, nverts):
		fsim[k] = pool.get(k)

	ind = np.argsort(fsim)
	fsim = np.take(fsim, ind, 0)
	sim = np.take(sim, ind, 0)

	iterations = 1
	simR = sim.copy() # reflection
	simE = sim.copy() # expansion
	simC = sim.copy() # contraction
	simFinal = sim.copy() # final
	fsimR = fsim.copy()
	fsimE = fsim.copy()
	fsimC = fsim.copy()
	fsimFinal = fsim.copy()


	while (fcalls[0] < maxfun and iterations < maxiter):
		if (np.max(np.ravel(np.abs(sim[1:] - sim[0]))) <= xtol and
			np.max(np.abs(fsim[0] - fsim[1:])) <= ftol):
			break

		xbar = np.add.reduce(sim[:-nactive], 0) / (nverts - nactive)

		state = np.zeros(nactive, dtype=int)
		actives = set(range(nactive))
		shrinks = []
		
		for i in range(nactive):
			# reflect bad points
			#xbar = (np.add.reduce(sim[:], 0)-sim[-i-1,:]) / N

			simR[-i - 1, :] = (1 + rho) * xbar - rho * sim[-i-1,:]
			applicator(i, simR[-i-1, :])
			state[i] = 1

		#print 'A'
		while len(actives)>0:

			#print 'ST', state
			i, curfsim = pool.get_any()
			#xbar = (np.add.reduce(sim[:], 0)-sim[-i-1,:]) / N

			if state[i] == 1: # point just got evaled after first reflection
				fsimR[-i-1] = curfsim
				if fsimR[-i-1] < fsim[0]: # expand
					simE[-i-1, :] = (1 + rho * chi) * xbar - rho * chi* sim[-i -1, :]
					applicator(i, simE[-i - 1, :])
					state[i] = 2
				else:
					if fsimR[-i-1] < fsim[-i-2]: # better than the next worst
						simFinal[-i-1, :] = simR[-i-1,:]
						fsimFinal[-i-1] = fsimR[-i-1]
						state[i]= -1 # we stopped after one reflection
						actives.remove(i)
					else: # contract
						if fsimR[-i-1] < fsim[-i-1]:
							simC[-i-1, :] = (1 + rho * psi) * xbar - rho * psi* sim[-i-1, :]
							applicator(i, simC[-i-1,:])
							state[i] = 3
						else:
							simC[-i- 1, :] = (1 - psi) * xbar + psi * sim[-i-1, :]
							applicator(i, simC[-i-1,:])
							state[i] = 4
			elif state[i] == 2: # the expansion step has been done 
				fsimE[-i-1] = curfsim
				#print 'TTT', fsimE[-i-1] , fsim[0]
				if fsimE[-i-1] < fsimR[-i-1]: # still better than the best
					simFinal[-i-1, :] = simE[-i-1,:]
					fsimFinal[-i-1] = fsimE[-i-1]
					state[i] = -2
					actives.remove(i)
				else:
					simFinal[-i-1, :] = simR[-i-1,:]
					fsimFinal[-i-1] = fsimR[-i-1]
					state[i] = -3
					actives.remove(i)
			elif state[i] == 3: # contraction1 just done
				fsimC[-i-1] = curfsim
				if fsimC[-i-1] < fsimR[-i-1]: # still better than the reflected
					simFinal[-i-1, :] = simC[-i-1,:]
					fsimFinal[-i-1] = fsimC[-i-1]
					state[i] = -4
					actives.remove(i)
				else:
					simFinal[-i-1, :] = sim[-i-1,:]
					fsimFinal[-i-1] = fsim[-i-1]
					shrinks.append(i)
					state[i] = -10
					actives.remove(i)
			elif state[i] == 4: # contraction2 just done
				fsimC[-i-1] = curfsim
				if fsimC[-i-1] < fsim[-i-1]: # still better than the reflected
					simFinal[-i-1, :] = simC[-i-1,:]
					fsimFinal[-i-1] = fsimC[-i-1]
					state[i] = -5
					actives.remove(i)
				else:
					simFinal[-i-1, :] = sim[-i-1,:]
					fsimFinal[-i-1] = fsim[-i-1]
					shrinks.append(i)
					state[i] = -10
					actives.remove(i)
			else:
				print 'Weird'
		for i in range(nactive):
			sim[-i-1,:] = simFinal[-i-1,:]
			fsim[-i-1] = fsimFinal[-i-1]
		#print 'ST1', state
		if len(shrinks) == nactive:
			print 'Shrink....'
			# no change was positive
			for j in one2np1:
				sim[j] = sim[0] + sigma * (sim[j] - sim[0])		
				applicator(j, sim[j, :])
			for j in one2np1:
				fsim[j] = pool.get(j)
		ind = np.argsort(fsim)
		sim = np.take(sim, ind, 0)
		fsim = np.take(fsim, ind, 0)

		if callback is not None:
			callback(sim[0])
		iterations += 1
		if retall:
			allvecs.append(sim)

	#allvecs = allcalls
	x = sim[0]
	fval = np.min(fsim)
	warnflag = 0

	if fcalls[0] >= maxfun:
		warnflag = 1
		msg = _status_message['maxfev']
		if disp:
			print('Warning: ' + msg)
	elif iterations >= maxiter:
		warnflag = 2
		msg = _status_message['maxiter']
		if disp:
			print('Warning: ' + msg)
	else:
		msg = _status_message['success']
		if disp:
			print(msg)
			print("		 Current function value: %f" % fval)
			print("		 Iterations: %d" % iterations)
			print("		 Function evaluations: %d" % fcalls[0])

	result = OptimizeResult(fun=fval, nit=iterations, nfev=fcalls[0],
					status=warnflag, success=(warnflag == 0), message=msg,
					x=x)
	if retall:
		result['allvecs'] = allvecs
	return result
				
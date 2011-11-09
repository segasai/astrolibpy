import glob, tables, numpy, Queue,psycopg2,threading,sys

def getConnection( db=None, driver=None, user=None,
						password=None, host=None):
	conn_str = "dbname=%s host=%s"%(db,host)
	if user is not None:
		conn_str = conn_str+ ' user=%s'%user
	if password is not None:
		conn_str = conn_str+ ' password=%s'%password
	conn = psycopg2.connect(conn_str)
	return conn					

def getCursor(conn, driver=None, preamb=None, notNamed=False):
	cur = conn.cursor()
	if preamb is not None:
		cur.execute(preamb)
	else:
		cur.execute('set cursor_tuple_fraction TO 1') 
		# this is required because otherwise PG may decide to execute a different plan
	if notNamed:
		return cur
	cur = conn.cursor(name='sqlutilcursor')
	cur.arraysize=100000
	return cur

def __inserter(desc, filename, qIn, endEvent):
	h5 = tables.openFile(filename, 'w')
	i=0;
	while(not endEvent.is_set()):
		try:
			tups = qIn.get(True,0.1)
			i+=1;
			#names='ra,dec,rmag,flag,psfr,aa,kk,airm,sky,camcol,typ'.split(',')
			#qOut.put(numpy.core.records.array(tups))
			names = [ x.name  for x in desc]
			if i==1:
				desc= {}
				desc1=[]
				for j,x, in enumerate(tups[0]):
					#desc[names[j]]=j.dtype
					desc1.append((names[j],type(x)))
				desc1=numpy.dtype(desc1)
				tab=h5.createTable('/','Table',desc1)
			tab.append(tups)
		except Queue.Empty:
			continue
		
	tab.close()
	h5.close()	

def get(query, filename, params=None, db="wsdb", driver="psycopg2", user=None,
						password=None, host='localhost', preamb=None,
						getConn=False, conn=None, maskNull=False):
	'''This program executes the sql query and saves the results
	in the HDF5 file (without storing all the results in memory)
	pg2hdf5.get('select ra,dec,d25 from rc3','rc3.h5')
	'''
	connSupplied = (conn is not None)
	if not connSupplied:
		conn = getConnection(db=db,driver=driver,user=user,password=password,
				host=host)
	try:
		cur = getCursor(conn, driver=driver, preamb=preamb)

		if params is None:
			res = cur.execute(query)
		else:
			res = cur.execute(query, params)

		qIn = Queue.Queue(1)
		endEvent = threading.Event()
		nrec = 0
		reslist = []
		try:
			while(True):
				tups = cur.fetchmany()
				if nrec==0:
					desc = cur.description
					proc = threading.Thread(target=__inserter, args = (desc, filename, qIn, endEvent))
					proc.start()

				if tups == []:
					break

				qIn.put(tups)
				nrec+=1
			endEvent.set()
		except BaseException:
			endEvent.set()
			proc.join(0.2) # notice that here the timeout is larger than the timeout
							# in the converter process
			if proc.is_alive():
				proc.terminate()
			raise
		proc.join()
		if reslist == []:
			nCols = len(cur.description)
			res = numpy.array([],
				dtype=numpy.dtype([('a%d'%i,'f') for i in range(nCols)])
								)				
		else:
			res = numpy.concatenate(reslist)
		res=[res[tmp] for tmp in res.dtype.names]
		if maskNull:
			for i in range(len(res)):
				if res[i].dtype==numpy.object:
					res[i]=res[i].astype(numpy.float)

	except BaseException:
		try:
			conn.rollback()
		except Exception:
			pass
		if not connSupplied:
			try:
				conn.close() # do not close if we were given the connection
			except:
				pass
		raise

	cur.close()
	conn.rollback()
	
	if not getConn:
		if not connSupplied:
			conn.close() # do not close if we were given the connection
		return res
	else:
		return conn,res


if __name__=='__main__':
	get(sys.argv[1],sys.argv[2],host='cappc118')
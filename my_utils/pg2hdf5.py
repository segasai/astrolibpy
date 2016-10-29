import glob, numpy, psycopg2,threading,sys
import h5py
try:
	from Queue import Queue, Empty
except:
	from queue import Queue, Empty
strLength = 20
__pgTypeHash = {16:bool,18:str,20:'i8',21:'i2',23:'i4',
		25: '|S%d' % strLength,
	700:'f4',701:'f8',
	1700:'f8',# numeric
	1114:'<M8[us]'  #timestamp
# 25:'|S%d'%strLength,
#	1042:'|S%d'%strLength,#character() 
#	1043:'|S%d'%strLength,#varchar
# I decided not supporting the strings
	} 

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
	cur.arraysize=1000000
	return cur

def __inserter(dtype, filename, qIn, endEvent):
	fp = h5py.File(filename, 'w')
	i=0;
	while(not endEvent.is_set()):
		try:
			tups = qIn.get(True,0.1)
			res=numpy.core.records.array(tups,dtype=dtype)
			if i==0:
				for _i, n in enumerate(res.dtype.names):
					fp.create_dataset(n,data=res[n],maxshape=(None,))
			else:
				newN = len(tups)
				N = len(fp[n])
				for _i, n in enumerate(res.dtype.names):
					fp[n].resize((N + newN,))
					fp[n][-newN:]=res[n]
			i+=1
		except Empty:
			continue
		
	fp.close()	

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

		qIn = Queue(1)
		endEvent = threading.Event()
		nrec = 0
		reslist = []
		try:
			while(True):
				tups = cur.fetchmany()
				if nrec==0:
					desc = cur.description
					typeCodes = [_tmp.type_code for _tmp in desc]
					colNames = [_tmp.name for _tmp in cur.description]
					dtype=numpy.dtype([(colNames[_i],__pgTypeHash[_t]) for _i,_t in enumerate(typeCodes)])				

					proc = threading.Thread(target=__inserter, args = (dtype, filename, qIn, endEvent))
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
	get(sys.argv[1],sys.argv[2],host='cappc127')
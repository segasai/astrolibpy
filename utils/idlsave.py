# -*- coding: utf-8 -*-
"""
To work with idlsave you need first to do 
from idlsave import idlsave

Then you can save and restore the variables 

idlsave.save('xx.sav','x,y',2,3)
exec (idlsave.restore('xx.sav'))

"""

import cPickle
import types


class idlsave:
	""" The class devoted for saving and restoring the python variable in
	files
	"""
	dhash={}
	def __init(self):
		pass
	
	@staticmethod
	def save(filename=None, names=None, *args):
		"""
		Saves your variable in a file, in such a way that you can easily
		retrieve them later. 
		Example:
		> x=2
		> y=[3,4,5]
		> idlsave.save("mydat.psav",'x,y',x,y)
		or 
		> exec(idlsave.save("myday.psav",'x,y'))
		Now you can leave python. You can later restore the x,y, variables 
		using idlsave.restore (see there for the doc)
		"""
		if len(args)==0:
			return "idlsave.save(\"%s\",\"%s\",%s)"%(filename,names,names)

		if type(names)==types.StringType:
			names=names.split(',')
		if len(names)!=len(args):
			raise Exception("The number of variable names should \
					be equal to the number of variables)")
		dhash={}
		for a in range(len(names)):
			dhash[names[a]]=args[a]
		f=open(filename,"w")
		cPickle.dump(dhash, f, 2)
		f.close()
		return None

	@staticmethod
	def restore(filename=None):
		"""Restores the variables stored in a file by idlsave.save routine
		Example: 
		> exec(idlsave.restore("mydat.psav"))
		
		Note that you MUST use this exact form exec(idlsave.restore(...))
		"""
		
		f=open(filename,"r")
		xx=cPickle.load(f)
		idlsave.dhash=xx
		f.close()
		buf=",".join(idlsave.dhash.iterkeys())
		if len(idlsave.dhash)==1:
			buf=buf+','
		buf=buf+"=idlsave.getallvars(filename=\"%s\")"%filename
		return buf

	@staticmethod
	def getallvars(filename=None):
		tup=tuple(a for a in idlsave.dhash.itervalues())
		dhash=None
		return tup

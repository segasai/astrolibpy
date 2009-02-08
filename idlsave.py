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
	dhash={}
	def __init(self):
		pass
	
	@staticmethod
	def save(filename=None, names=None, *args):

		if type(names)==types.StringType:
			names=names.split(',')
		if len(names)!=len(args):
			raise Exception("The number of variable names should \
					be equal to the number of variables)")
		dhash={}
		for a in range(len(names)):
			dhash[names[a]]=args[a]
		f=open(filename,"w")
		cPickle.dump(dhash,f)
		f.close()
		return None

	@staticmethod
	def restore(filename=None):
		f=open(filename,"r")
		xx=cPickle.load(f)
		idlsave.dhash=xx
		f.close()
		buf=""
		for a in idlsave.dhash.iterkeys():
			buf=buf+a+","
		buf=buf[0:-1]
		buf=buf+"=idlsave.getallvars(filename=\"%s\")"%filename
		return buf

	@staticmethod
	def getallvars(filename=None):
		tup=tuple([None for i in range(len(idlsave.dhash)) ] )
		i=0
		tup=tuple(a for a in idlsave.dhash.itervalues())
		dhash=None
		return tup

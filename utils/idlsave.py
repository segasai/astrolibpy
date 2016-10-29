# -*- coding: utf-8 -*-
# Copyright (C) 2009-2010 Sergey Koposov
# This file is part of astrolibpy
#
#	 astrolibpy is free software: you can redistribute it and/or modify
#	 it under the terms of the GNU General Public License as published by
#	 the Free Software Foundation, either version 3 of the License, or
#	 (at your option) any later version.
#
#	astrolibpy is distributed in the hope that it will be useful,
#	 but WITHOUT ANY WARRANTY; without even the implied warranty of
#	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	 GNU General Public License for more details.
#
#	 You should have received a copy of the GNU General Public License
#	 along with astrolibpy.	 If not, see <http://www.gnu.org/licenses/>.


"""
To work with idlsave you need first to do 
from idlsave import idlsave

Then you can save and restore the variables: 
> x=2
> y=3
> idlsave.save('xx.sav','x,y', x, y)
OR
> exec(idlsave.save('xx.sav','x,y'))

To restore you need to do:
> exec (idlsave.restore('xx.sav'))

"""
from __future__ import print_function
try:
	import cPickle as pickle
except ImportError:
	import pickle
import types
import re
import struct


class idlsave:
	""" The class designed for saving and restoring python variables to
	files
	"""
	dhash = {}

	def __init(self):
		pass

	@staticmethod
	def versionId(version):
		return 'IDLSAVE-v%04d' % version

	@staticmethod
	def parseVersion(s):
		m = re.match('IDLSAVE-v(\d\d\d\d)', s)
		if m is None:
			return 1
		else:
			return int(m.group(1))

	@staticmethod
	def __cureString(s):
		return s.replace('\n', '').replace(' ', '').replace('\t', '')

	@staticmethod
	def __splitString(s):
		return idlsave.__cureString(s).split(',')

	@staticmethod
	def save(filename=None, names=None, *args, **kw):
		"""
		Saves your variable in a file, in such a way that you can easily
		retrieve them later. 
		Example:
		> x = 2
		> y = [3,4,5]
		> idlsave.save("mydat.psav", 'x,y', x, y)
		or 
		> exec(idlsave.save("myday.psav",'x,y'))
		Now you can leave python. You can later restore the x,y, variables 
		using idlsave.restore
		
		Storage format for version 2 files is the following:
			The first 8 bytes store the offset of the pickled dictionary of 
			offsets. The dictionary of offsets stores the
			offset of each pickled variable in the file.
			The keys in the dictionary are variable names.
		"""

		if len(args) == 0:
			return "idlsave.save(\"%s\",\"%s\",%s)" % (filename, idlsave.__cureString(names), names)

		if isinstance(names, str):
			names = idlsave.__splitString(names)
		if len(names) != len(args):
			raise Exception("The number of variable names should \
					be equal to the number of variables)")
		f = open(filename, "wb")
		version = kw.get('version', 2)
		curhash = {}
		for a in range(len(names)):
			curhash[names[a]] = args[a]

		if version == 1:
			pickle.dump(curhash, f, pickle.HIGHEST_PROTOCOL)
		elif version == 2:
			f.write(idlsave.versionId(version).encode('ascii'))
			headlen1 = f.tell()
			f.write(struct.pack('!q', 0))
			offsets = dict([(name, 0) for name in names])
			for name in names:
				offsets[name] = f.tell()
				pickle.dump(curhash[name], f, pickle.HIGHEST_PROTOCOL)
			offOffs = f.tell()
			pickle.dump(offsets, f, pickle.HIGHEST_PROTOCOL)
			f.seek(headlen1)
			f.write(struct.pack('!q', offOffs))
		f.close()
		del curhash
		return None

	@staticmethod
	def restore(filename=None, names=None, asdict=False, version=None, printVars=False):
		"""Restores the variables stored in a file by idlsave.save routine
		Example: 
		> exec(idlsave.restore("mydat.psav"))
		Note that you MUST use this exact form exec(idlsave.restore(...))
		Or you can restore only the variables that you want:
		> ra, dec = idlsave.restore("mysav.psav", "ra,dec")
		"""
		f = open(filename, "rb")
		vid = idlsave.versionId(1)
		try:
			prefix = f.read(len(vid)).decode('ascii')
		except:	
			version = 1
		if version is None:
			version = idlsave.parseVersion(prefix)
		if version == 1:
			f.seek(0)  # there is no header in the version1 files

		if version == 1:
			xx = pickle.load(f)
			idlsave.dhash = xx
			f.close()
			if names is None:
				if asdict:
					ret = idlsave.dhash.copy()
					del idlsave.dhash
					return ret
				buf = ",".join(idlsave.dhash.iterkeys())
				if len(idlsave.dhash) == 1:
					buf = buf + ','
				buf = buf + \
					"=idlsave.getallvars(printVars=%s)" % (str(printVars))
				return buf
			else:
				names = idlsave.__splitString(names)
				res = [idlsave.dhash[a] for a in names]
				del idlsave.dhash
				return res
		elif version == 2:
			offOff = struct.unpack('!q', f.read(8))[0]
			f.seek(offOff)
			offsets = pickle.load(f)
			if names is None:
				names1 = offsets.keys()
			else:
				names1 = idlsave.__splitString(names)
			hash = {}
			for name in names1:
				off = offsets[name]
				try:
					f.seek(off)
					hash[name] = pickle.load(f)
				except UnicodeDecodeError:
					f.seek(off)
					hash[name] = pickle.load(f, encoding='latin1')
			if asdict:
				return hash
			idlsave.dhash = hash
			f.close()
			if names is None:
				buf = ",".join(idlsave.dhash.keys())
				if len(idlsave.dhash) == 1:
					buf = buf + ','
				buf = buf + "=idlsave.getallvars(%s)" % str(printVars)
				return buf	# return the string for exec
			else:
				res = [idlsave.dhash[a] for a in names1]
				del idlsave.dhash
				return res

	@staticmethod
	def getallvars(printVars=False):
		tup = tuple(a for a in idlsave.dhash.values())
		if printVars:
			print (','.join([k for k in idlsave.dhash.keys()]))
		del idlsave.dhash
		return tup

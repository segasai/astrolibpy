# Copyright (C) 2009-2010 Sergey Koposov
# This file is part of astrolibpy
#
#    astrolibpy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#   astrolibpy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with astrolibpy.  If not, see <http://www.gnu.org/licenses/>.


import scipy.spatial.kdtree, numpy
from scipy import version
from numpy import sin,cos,deg2rad,rad2deg,arcsin

scipy_version  = float('.'.join(version.version.split('.')[0:2]))

def match_lists(ra1, dec1, ra2, dec2, dist):
	# crossmatches the list of objects (ra1,dec1) with
	# another list of objects (ra2,dec2) with the dist matching radius
	# the routine returns the distance to the neighbor and the list 
	# of indices of the neighbor. Everything is in degrees.
	# if no match is found the distance is NaN.
	# Example 
	# dist,ind=match(ra1,dec1,ra2,dec2)
	# goodmatch_ind = numpy.isfinite(dist)
	# plot(ra1[goodmatch_ind],ra2[ind][goodmatch_ind])
	cosd = lambda x : cos(deg2rad(x))
	sind = lambda x : sin(deg2rad(x))
	mindist = 2 * sind(dist/2)	
	getxyz = lambda r, d: [cosd(r)*cosd(d), sind(r)*cosd(d), sind(d)]
	xyz1 = numpy.array(getxyz(ra1, dec1))
	xyz2 = numpy.array(getxyz(ra2, dec2))
	
	if scipy_version<0.8:
	# If old scipy version is detected then we use KDTree instead of 
	# cKDTtree because there is a bug in the cKDTree
	# http://projects.scipy.org/scipy/ticket/1178
		tree2 = scipy.spatial.KDTree(xyz2.T)
	else:
		tree2 = scipy.spatial.cKDTree(xyz2.T)
	ret = tree2.query(xyz1.T, 1, 0, 2, mindist)
	dist, ind = ret
	dist = rad2deg(2*arcsin(dist/2))
	return dist, ind
	
	
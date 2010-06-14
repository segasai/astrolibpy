import scipy.spatial.kdtree, numpy
from numpy import sin,cos,deg2rad,rad2deg,arcsin


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

# At the moment I'm using Python version of the KDTree instead of 
# cKDTtree because there is a bug in the cKDTree
# http://projects.scipy.org/scipy/ticket/1178
	tree2 = scipy.spatial.KDTree(xyz2.transpose())
	ret = tree2.query(xyz1.transpose(), 1, 0, 2, mindist)
	dist, ind = ret
	dist = rad2deg(2*arcsin(dist/2))
	return dist, ind
	
	
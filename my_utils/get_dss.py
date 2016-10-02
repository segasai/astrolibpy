from __future__ import print_function
import aplpy
try:
	from urllib2 import urlopen
except:
	from urllib.request import urlopen
try:
	from StringIO import StringIO
except:
	from io import BytesIO as StringIO
		
import astropy.io.fits as pyfits
import matplotlib.pyplot  as plt
def get_dss(ra, dec, survey='all', radius = 15, debug=False, noerase=False,
		subplot=(1,1,1)):
	"""
	radius is in arcmins
	survey: 'all', 'poss2ukstu', 'dss1'.. 
	"""
	#survey='poss2ukstu'
	url='http://archive.stsci.edu/cgi-bin/dss_search?v=%s&r=%fd&d=%fd&e=J2000&h=%f&w=%f&f=fits&c=none&fov=NONE&v3='%(survey,ra,dec,radius,radius)
	
	response = urlopen(url)
	if debug:
		print( url)
	html = response.read()
	f=StringIO(html)
	
	dat=pyfits.open(f)
	if not noerase:
		plt.clf()
	gc=aplpy.FITSFigure(dat,figure=plt.gcf(),subplot=subplot)
	gc.set_tick_labels_format('ddd.dddd','ddd.dddd')
	gc.show_grayscale()
	return gc

	
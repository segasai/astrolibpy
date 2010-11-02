import urllib2,StringIO,aplpy,pyfits
import matplotlib.pyplot  as plt
def get_dss(ra, dec, survey='all', radius = 15, debug=False):
	"""
	radius is in arcmins
	survey: 'all', 'poss2ukstu', 'dss1'.. 
	"""
	#survey='poss2ukstu'
	url='http://archive.stsci.edu/cgi-bin/dss_search?v=%s&r=%fd&d=%fd&e=J2000&h=%f&w=%f&f=fits&c=none&fov=NONE&v3='%(survey,ra,dec,radius,radius)
	
	response = urllib2.urlopen(url)
	if debug:
		print url
	html = response.read()
	f=StringIO.StringIO(html)
	
	dat=pyfits.open(f)
	plt.clf()
	gc=aplpy.FITSFigure(dat,figure=plt.gcf())
	gc.set_tick_labels_format('ddd.dddd','ddd.dddd')
	gc.show_grayscale()
	

	
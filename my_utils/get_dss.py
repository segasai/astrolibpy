import urllib2,StringIO,aplpy,pyfits
import matplotlib.pyplot  as plt
def get_dss(ra,dec):
	url='http://archive.stsci.edu/cgi-bin/dss_search?v=poss2ukstu_red&r=%fd&d=%fd&e=J2000&h=15.0&w=15.0&f=fits&c=none&fov=NONE&v3='%(ra,dec)
	
	response = urllib2.urlopen(url)
	html = response.read()
	f=StringIO.StringIO(html)
	
	dat=pyfits.open(f)
	plt.clf()
	gc=aplpy.FITSFigure(dat,figure=plt.gcf())
	gc.set_tick_labels_format('ddd.dddd','ddd.dddd')
	gc.show_grayscale()
	

	
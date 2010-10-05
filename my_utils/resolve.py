import sys
from suds.client import Client
from lxml import etree

def resolve(objectName):
	url = 'http://cdsws.u-strasbg.fr/axis/services/Sesame?wsdl'
	client = Client(url)    
	xml=client.service.SesameXML(objectName)           
	tree = etree.fromstring(xml.encode('utf-8'))
	# take the first resolver
	pathRa = tree.xpath('/Sesame/Target/Resolver[1]/jradeg')
	pathDec = tree.xpath('/Sesame/Target/Resolver[1]/jdedeg')
	if len(pathRa)==0:
		return []
	ra=float(pathRa[0].text)
	dec=float(pathDec[0].text)
	return ra,dec


if __name__=='__main__':
	res = resolve(sys.argv[1])
	if len(res)>0:
		print res[0],res[1]
	else:
		print 'Not found'

import sys
#from suds.client import Client
from zeep import Client
from lxml import etree

def resolve(objectName):
    """    R
esolve the object by name using CDS
    Example:
    >> resolve.resolve('M31')
    (10.684708329999999, 41.268749999999997)

    Requires the following modules:
        suds, lxml
    """
    url = 'http://cdsws.u-strasbg.fr/axis/services/Sesame?wsdl'
    client = Client(url)    
    xml=client.service.SesameXML(objectName)           
    tree = etree.fromstring(xml.encode('utf-8'))
    # take the first resolver
    success=True
    for i in range(4):
        pathRa = tree.xpath('/Sesame/Target/Resolver[%d]/jradeg'%i)
        pathDec = tree.xpath('/Sesame/Target/Resolver[%d]/jdedeg'%i)
        if len(pathRa)!=0:
            success=True
            break
    if not success:
        return []
    ra=float(pathRa[0].text)
    dec=float(pathDec[0].text)
    return ra,dec


if __name__=='__main__':
    res = resolve(sys.argv[1])
    if len(res)>0:
        print (res[0],res[1])
    else:
        print ('Not found')

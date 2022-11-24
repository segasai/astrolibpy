import sys
from urllib.parse import quote_plus
from lxml import etree
from urllib.request import urlopen


def resolve(objectName):
    """    R
esolve the object by name using CDS
    Example:
    >> resolve.resolve('M31')
    (10.684708329999999, 41.268749999999997)

    Requires the following modules:
        suds, lxml
    """
    url = 'http://vizier.u-strasbg.fr/cgi-bin/Sesame/-ox/?%s' % quote_plus(
        objectName)
    f = urlopen(url)
    result = f.read()
    tree = etree.fromstring(result)
    # take the first resolver
    success = True
    for i in range(4):
        pathRa = tree.xpath('/Sesame/Target/Resolver[%d]/jradeg' % i)
        pathDec = tree.xpath('/Sesame/Target/Resolver[%d]/jdedeg' % i)
        if len(pathRa) != 0:
            success = True
            break
    if not success:
        return []
    ra = float(pathRa[0].text)
    dec = float(pathDec[0].text)
    return ra, dec


if __name__ == '__main__':
    res = resolve(sys.argv[1])
    if len(res) > 0:
        print(res[0], res[1])
    else:
        print('Not found')

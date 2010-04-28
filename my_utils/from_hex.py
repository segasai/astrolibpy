import numpy, re
def from_hex(arr):
	r=re.compile('(\-?)(.+):(.+):(.+)')
	ret=[]
	for a in arr:
		m = r.search(a)

		sign = m.group(1)=='-'
		if sign:
			sign=-1
		else:
			sign=1
		i1 = int(m.group(2))
		i2 = int(m.group(3))
		i3 = float(m.group(4))
		val = sign*(int(i1)+int(i2)/60.+(float(i3))/3600.)
		ret.append(val)
	return numpy.array(ret)
		
		
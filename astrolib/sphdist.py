from numpy import deg2rad,rad2deg,sin,cos,sqrt,arcsin

def sphdist (ra1, dec1, ra2, dec2):
	"""measures the spherical distance in degrees
		The input has to be in degrees too
	"""
	dec1_r = deg2rad(dec1)
	dec2_r = deg2rad(dec2)
	return 2 *\
		rad2deg \
			(
			arcsin
				(
				sqrt
					(
						(
							sin((dec1_r - dec2_r) / 2)
						)**2
						+  
						cos(dec1_r) * cos(dec2_r) *
						(
							sin((deg2rad(ra1 - ra2)) / 2)
						)**2
					)
				)
			)
	
def sphdist_fast(ra1,dec1,ra2,dec2):
	import numexpr
	return numexpr.evaluate('2*57.295779513082323*(arcsin(sqrt((sin(0.017453292519943295*(dec1-dec2)/2))**2+cos(0.017453292519943295*dec1)*cos(0.017453292519943295*dec2)*(sin(0.017453292519943295*((ra1-ra2))/2))**2)))')
	
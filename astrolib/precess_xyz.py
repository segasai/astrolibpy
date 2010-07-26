from numpy import sin, cos, arctan2, sqrt, arcsin
from precess import precess

def precess_xyz(x, y, z, equinox1, equinox2):
#+
# NAME:
#	PRECESS_XYZ
#
# PURPOSE:
#	Precess equatorial geocentric rectangular coordinates.
#
# CALLING SEQUENCE:
#	precess_xyz, x, y, z, equinox1, equinox2
#
# INPUT/OUTPUT:
#	x,y,z: scalars or vectors giving heliocentric rectangular coordinates
#              THESE ARE CHANGED UPON RETURNING.
# INPUT:
#	EQUINOX1: equinox of input coordinates, numeric scalar
#       EQUINOX2: equinox of output coordinates, numeric scalar
#
# OUTPUT:
#	x,y,z are changed upon return
#
# NOTES:
#   The equatorial geocentric rectangular coords are converted
#      to RA and Dec, precessed in the normal way, then changed
#      back to x, y and z using unit vectors.
#
#EXAMPLE:
#	Precess 1950 equinox coords x, y and z to 2000.
#	IDL> precess_xyz,x,y,z, 1950, 2000
#
#HISTORY:
#	Written by P. Plait/ACC March 24 1999
#	   (unit vectors provided by D. Lindler)
#       Use /Radian call to PRECESS     W. Landsman     November 2000
#       Use two parameter call to ATAN   W. Landsman    June 2001
#-
#check inputs
   
   #take input coords and convert to ra and dec (in radians)
   
   ra = arctan2(y, x)
   _del = sqrt(x * x + y * y + z * z)  #magnitude of distance to Sun
   dec = arcsin(z / _del)
   
   #   precess the ra and dec
   ra,dec = precess(ra, dec, equinox1, equinox2, radian=True)
   
   #convert back to x, y, z
   xunit = cos(ra) * cos(dec)
   yunit = sin(ra) * cos(dec)
   zunit = sin(dec)
   
   x = xunit * _del
   y = yunit * _del
   z = zunit * _del
   
   return x,y,z



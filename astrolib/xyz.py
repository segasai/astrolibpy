from numpy import sqrt, pi, cos, sin
from precess_xyz import precess_xyz

def xyz(date, equinox=None):
   """
    NAME:
          XYZ
    PURPOSE:
          Calculate geocentric X,Y, and Z  and velocity coordinates of the Sun
    EXPLANATION:
          Calculates geocentric X,Y, and Z vectors and velocity coordinates
          (dx, dy and dz) of the Sun.   (The positive X axis is directed towards
          the equinox, the y-axis, towards the point on the equator at right
          ascension 6h, and the z axis toward the north pole of the equator).
          Typical position accuracy is <1e-4 AU (15000 km).
   
    CALLING SEQUENCE:
          XYZ, date, x, y, z, [ xvel, yvel, zvel, EQUINOX = ]
   
    INPUT:
          date: reduced julian date (=JD - 2400000), scalar or vector
   
    OUTPUT:
          x,y,z: scalars or vectors giving heliocentric rectangular coordinates
                    (in A.U) for each date supplied.    Note that sqrt(x^2 + y^2
                    + z^2) gives the Earth-Sun distance for the given date.
          xvel, yvel, zvel: velocity vectors corresponding to X, Y and Z.
   
    OPTIONAL KEYWORD INPUT:
          EQUINOX: equinox of output. Default is 1950.
   
    EXAMPLE:
          What were the rectangular coordinates and velocities of the Sun on
          Jan 22, 1999 0h UT (= JD 2451200.5) in J2000 coords? NOTE:
          Astronomical Almanac (AA) is in TDT, so add 64 seconds to
          UT to convert.
   
          IDL> xyz,51200.5+64.d/86400.d,x,y,z,xv,yv,zv,equinox = 2000
   
          Compare to Astronomical Almanac (1999 page C20)
                      X  (AU)        Y  (AU)     Z (AU)
          XYZ:      0.51456871   -0.76963263  -0.33376880
          AA:       0.51453130   -0.7697110   -0.3337152
          abs(err): 0.00003739    0.00007839   0.00005360
          abs(err)
              (km):   5609          11759         8040
   
          NOTE: Velocities in AA are for Earth/Moon barycenter
                (a very minor offset) see AA 1999 page E3
                     X VEL (AU/DAY) YVEL (AU/DAY)   Z VEL (AU/DAY)
          XYZ:      -0.014947268   -0.0083148382    -0.0036068577
          AA:       -0.01494574    -0.00831185      -0.00360365
          abs(err):  0.000001583    0.0000029886     0.0000032077
          abs(err)
           (km/sec): 0.00265        0.00519          0.00557
   
    PROCEDURE CALLS:
          PRECESS_XYZ
    REVISION HISTORY
          Original algorithm from Almanac for Computers, Doggett et al. USNO 1978
          Adapted from the book Astronomical Photometry by A. Henden
          Written  W. Landsman   STX       June 1989
          Correct error in X coefficient   W. Landsman HSTX  January 1995
          Added velocities, more terms to positions and EQUINOX keyword,
             some minor adjustments to calculations
             P. Plait/ACC March 24, 1999
   """

   picon = pi / 180.0e0
   t = (date - 15020.0e0) / 36525.0e0         #Relative Julian century from 1900
   
   # NOTE: longitude arguments below are given in *equinox* of date.
   #   Precess these to equinox 1950 to give everything an even footing.
   #   Compute argument of precession from equinox of date back to 1950
   pp = (1.396041e0 + 0.000308e0 * (t + 0.5e0)) * (t - 0.499998e0)
   
   # Compute mean solar longitude, precessed back to 1950
   el = 279.696678e0 + 36000.76892e0 * t + 0.000303e0 * t * t - pp
   
   # Compute Mean longitude of the Moon
   c = 270.434164e0 + 480960.e0 * t + 307.883142e0 * t - 0.001133e0 * t * t - pp
   
   # Compute longitude of Moon's ascending node
   n = 259.183275e0 - 1800.e0 * t - 134.142008e0 * t + 0.002078e0 * t * t - pp
   
   # Compute mean solar anomaly
   g = 358.475833e0 + 35999.04975e0 * t - 0.00015e0 * t * t
   
   # Compute the mean jupiter anomaly
   j = 225.444651e0 + 2880.0e0 * t + 154.906654e0 * t * t
   
   # Compute mean anomaly of Venus
   v = 212.603219e0 + 58320.e0 * t + 197.803875e0 * t + 0.001286e0 * t * t
   
   # Compute mean anomaly of Mars
   m = 319.529425e0 + 19080.e0 * t + 59.8585e0 * t + 0.000181e0 * t * t
   
   # Convert degrees to radians for trig functions
   el = el * picon
   g = g * picon
   j = j * picon
   c = c * picon
   v = v * picon
   n = n * picon
   m = m * picon
   
   # Calculate X,Y,Z using trigonometric series
   x = 0.999860e0 * cos(el) - 0.025127e0 * cos(g - el) + 0.008374e0 * cos(g + el) + 0.000105e0 * cos(g + g + el) + 0.000063e0 * t * cos(g - el) + 0.000035e0 * cos(g + g - el) - 0.000026e0 * sin(g - el - j) - 0.000021e0 * t * cos(g + el) + 0.000018e0 * sin(2.e0 * g + el - 2.e0 * v) + 0.000017e0 * cos(c) - 0.000014e0 * cos(c - 2.e0 * el) + 0.000012e0 * cos(4.e0 * g + el - 8.e0 * m + 3.e0 * j) - 0.000012e0 * cos(4.e0 * g - el - 8.e0 * m + 3.e0 * j) - 0.000012e0 * cos(g + el - v) + 0.000011e0 * cos(2.e0 * g + el - 2.e0 * v) + 0.000011e0 * cos(2.e0 * g - el - 2.e0 * j)
   
   
   y = 0.917308e0 * sin(el) + 0.023053e0 * sin(g - el) + 0.007683e0 * sin(g + el) + 0.000097e0 * sin(g + g + el) - 0.000057e0 * t * sin(g - el) - 0.000032e0 * sin(g + g - el) - 0.000024e0 * cos(g - el - j) - 0.000019e0 * t * sin(g + el) - 0.000017e0 * cos(2.e0 * g + el - 2.e0 * v) + 0.000016e0 * sin(c) + 0.000013e0 * sin(c - 2.e0 * el) + 0.000011e0 * sin(4.e0 * g + el - 8.e0 * m + 3.e0 * j) + 0.000011e0 * sin(4.e0 * g - el - 8.e0 * m + 3.e0 * j) - 0.000011e0 * sin(g + el - v) + 0.000010e0 * sin(2.e0 * g + el - 2.e0 * v) - 0.000010e0 * sin(2.e0 * g - el - 2.e0 * j)
   
   
   z = 0.397825e0 * sin(el) + 0.009998e0 * sin(g - el) + 0.003332e0 * sin(g + el) + 0.000042e0 * sin(g + g + el) - 0.000025e0 * t * sin(g - el) - 0.000014e0 * sin(g + g - el) - 0.000010e0 * cos(g - el - j)
   
   #Precess_to new equator?
   if equinox is not None:   
      x, y, z = precess_xyz(x, y, z, 1950, equinox)
   
   xvel = -0.017200e0 * sin(el) - 0.000288e0 * sin(g + el) - 0.000005e0 * sin(2.e0 * g + el) - 0.000004e0 * sin(c) + 0.000003e0 * sin(c - 2.e0 * el) + 0.000001e0 * t * sin(g + el) - 0.000001e0 * sin(2.e0 * g - el)
   
   yvel = 0.015780 * cos(el) + 0.000264 * cos(g + el) + 0.000005 * cos(2.e0 * g + el) + 0.000004 * cos(c) + 0.000003 * cos(c - 2.e0 * el) - 0.000001 * t * cos(g + el)
   
   zvel = 0.006843 * cos(el) + 0.000115 * cos(g + el) + 0.000002 * cos(2.e0 * g + el) + 0.000002 * cos(c) + 0.000001 * cos(c - 2.e0 * el)
   
   #Precess to new equator?
   
   if equinox is not None:   
      xvel, yvel, zvel = precess_xyz(xvel, yvel, zvel, 1950, equinox)
   
   return x, y, z, xvel, yvel, zvel


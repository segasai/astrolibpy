# -*- coding: utf-8 -*-
from numpy import array, sin, cos, zeros, pi, zeros

def premat(equinox1, equinox2, fk4=False):
   """
    NAME:
          PREMAT
    PURPOSE:
          Return the precession matrix needed to go from EQUINOX1 to EQUINOX2.
    EXPLANTION:
          This matrix is used by the procedures PRECESS and BARYVEL to precess
          astronomical coordinates
   
    CALLING SEQUENCE:
          matrix = PREMAT( equinox1, equinox2, [ /FK4 ] )
   
    INPUTS:
          EQUINOX1 - Original equinox of coordinates, numeric scalar.
          EQUINOX2 - Equinox of precessed coordinates.
   
    OUTPUT:
         matrix - double precision 3 x 3 precession matrix, used to precess
                  equatorial rectangular coordinates
   
    OPTIONAL INPUT KEYWORDS:
          /FK4   - If this keyword is set, the FK4 (B1950.0) system precession
                  angles are used to compute the precession matrix.   The
                  default is to use FK5 (J2000.0) precession angles
   
    EXAMPLES:
          Return the precession matrix from 1950.0 to 1975.0 in the FK4 system
   
          IDL> matrix = PREMAT( 1950.0, 1975.0, /FK4)
   
    PROCEDURE:
          FK4 constants from "Computational Spherical Astronomy" by Taff (1983),
          p. 24. (FK4). FK5 constants from "Astronomical Almanac Explanatory
          Supplement 1992, page 104 Table 3.211.1.
   
    REVISION HISTORY
          Written, Wayne Landsman, HSTX Corporation, June 1994
          Converted to IDL V5.0   W. Landsman   September 1997
   """

   deg_to_rad = pi / 180.0e0
   sec_to_rad = deg_to_rad / 3600.e0
   
   t = 0.001e0 * (equinox2 - equinox1)
   
   if not fk4:   
      st = 0.001e0 * (equinox1 - 2000.e0)
      #  Compute 3 rotation angles
      a = sec_to_rad * t * (23062.181e0 + st * (139.656e0 + 0.0139e0 * st) + t * (30.188e0 - 0.344e0 * st + 17.998e0 * t))
      
      b = sec_to_rad * t * t * (79.280e0 + 0.410e0 * st + 0.205e0 * t) + a
      
      c = sec_to_rad * t * (20043.109e0 - st * (85.33e0 + 0.217e0 * st) + t * (-42.665e0 - 0.217e0 * st - 41.833e0 * t))
      
   else:   
      
      st = 0.001e0 * (equinox1 - 1900.e0)
      #  Compute 3 rotation angles
      
      a = sec_to_rad * t * (23042.53e0 + st * (139.75e0 + 0.06e0 * st) + t * (30.23e0 - 0.27e0 * st + 18.0e0 * t))
      
      b = sec_to_rad * t * t * (79.27e0 + 0.66e0 * st + 0.32e0 * t) + a
      
      c = sec_to_rad * t * (20046.85e0 - st * (85.33e0 + 0.37e0 * st) + t * (-42.67e0 - 0.37e0 * st - 41.8e0 * t))
      
   
   sina = sin(a)
   sinb = sin(b)
   sinc = sin(c)
   cosa = cos(a)
   cosb = cos(b)
   cosc = cos(c)
   
   r = zeros((3, 3))
   r[0,:] = array([cosa * cosb * cosc - sina * sinb, sina * cosb + cosa * sinb * cosc, cosa * sinc])
   r[1,:] = array([-cosa * sinb - sina * cosb * cosc, cosa * cosb - sina * sinb * cosc, -sina * sinc])
   r[2,:] = array([-cosb * sinc, -sinb * sinc, cosc])
   
   return r


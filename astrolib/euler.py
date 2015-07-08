from numpy import array, sin, cos, pi, deg2rad, rad2deg, arctan2, arcsin, minimum

def euler(ai, bi, select=1, fk4=False):
   """
    NAME:
        EULER
    PURPOSE:
        Transform between Galactic, celestial, and ecliptic coordinates.
    EXPLANATION:
        Use the procedure ASTRO to use this routine interactively
   
    CALLING SEQUENCE:
         AO, BO = EULER(AI, BI, [SELECT=1, FK4=False])
   
    INPUTS:
          AI - Input Longitude in DEGREES, scalar or vector.  If only two
                  parameters are supplied, then  AI and BI will be modified to
                  contain the output longitude and latitude.
          BI - Input Latitude in DEGREES
   
    OPTIONAL INPUT:
          SELECT - Integer (1-6) specifying type of coordinate transformation.
   
         SELECT   From          To        |   SELECT      From            To
          1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec
          2     Galactic       RA-DEC     |     5       Ecliptic      Galactic
          3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic
   
         If not supplied as a parameter or keyword, then EULER will prompt for
         the value of SELECT
         Celestial coordinates (RA, Dec) should be given in equinox J2000
         unless the /FK4 keyword is set.
    OUTPUTS:
          AO - Output Longitude in DEGREES
          BO - Output Latitude in DEGREES
   
    INPUT KEYWORD:
          /FK4 - If this keyword is set and non-zero, then input and output
                celestial and ecliptic coordinates should be given in equinox
                B1950.
          /SELECT  - The coordinate conversion integer (1-6) may alternatively be
                 specified as a keyword
    NOTES:
          EULER was changed in December 1998 to use J2000 coordinates as the
          default, ** and may be incompatible with earlier versions***.
    REVISION HISTORY:
          Written W. Landsman,  February 1987
          Adapted from Fortran by Daryl Yentis NRL
          Converted to IDL V5.0   W. Landsman   September 1997
          Made J2000 the default, added /FK4 keyword  W. Landsman December 1998
          Add option to specify SELECT as a keyword W. Landsman March 2003
   """

   #   J2000 coordinate conversions are based on the following constants
   #   (see the Hipparcos explanatory supplement).
   #  eps = 23.4392911111d              Obliquity of the ecliptic
   #  alphaG = 192.85948d               Right Ascension of Galactic North Pole
   #  deltaG = 27.12825d                Declination of Galactic North Pole
   #  lomega = 32.93192d                Galactic longitude of celestial equator
   #  alphaE = 180.02322d              Ecliptic longitude of Galactic North Pole
   #  deltaE = 29.811438523d            Ecliptic latitude of Galactic North Pole
   #  Eomega  = 6.3839743d              Galactic longitude of ecliptic equator
   
   if fk4:   
      equinox = '(B1950)'
      psi = array ([0.57595865315e0, 4.9261918136e0, 0.00000000000e0, 0.0000000000e0, 0.11129056012e0, 4.7005372834e0])
      stheta = array ([0.88781538514e0, -0.88781538514e0, 0.39788119938e0, -0.39788119938e0, 0.86766174755e0, -0.86766174755e0])
      ctheta = array([0.46019978478e0, 0.46019978478e0, 0.91743694670e0, 0.91743694670e0, 0.49715499774e0, 0.49715499774e0])
      phi = array([4.9261918136e0, 0.57595865315e0, 0.0000000000e0, 0.00000000000e0, 4.7005372834e0, 0.11129056012e0])
   else:   
      equinox = '(J2000)'
      psi = array([0.57477043300e0, 4.9368292465e0, 0.00000000000e0, 0.0000000000e0, 0.11142137093e0, 4.71279419371e0])
      stheta = array([0.88998808748e0, -0.88998808748e0, 0.39777715593e0, -0.39777715593e0, 0.86766622025e0, -0.86766622025e0])
      ctheta = array([0.45598377618e0, 0.45598377618e0, 0.91748206207e0, 0.91748206207e0, 0.49714719172e0, 0.49714719172e0])
      phi = array([4.9368292465e0, 0.57477043300e0, 0.0000000000e0, 0.00000000000e0, 4.71279419371e0, 0.11142137093e0])
   if select not in [1,2,3,4,5,6]:
      raise ValueError('Select parameter should be an integer between 1 and 6')
   i = select - 1
   b = deg2rad(bi)
   cb = cos(b)
   sb = sin(b)
   del b
   a = deg2rad(ai) - phi[i]
   cbsa = cb * sin(a)
   b = -stheta[i] * cbsa + ctheta[i] * sb
   bo = rad2deg(arcsin(minimum(b, 1.0)))
   del b
   a = arctan2(ctheta[i] * cbsa + stheta[i] * sb, cb * cos(a))
   del cb, cbsa, sb
   ao = rad2deg(((a + psi[i] + 4 * pi ) % (2 * pi) ) )

   return (ao, bo)

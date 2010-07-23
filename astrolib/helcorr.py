from numpy import *
from baryvel import baryvel
from daycnv import daycnv
from precess import precess
from helio_jd import helio_jd
_radeg = 180.0 / pi

def helcorr(obs_long, obs_lat, obs_alt, ra2000, dec2000, jd, debug=False):

#calculates heliocentric Julian date, baricentric and heliocentric radial
#velocity corrections from:
#
#INPUT:
#<OBSLON> Longitude of observatory (degrees, western direction is positive)
#<OBSLAT> Latitude of observatory (degrees)
#<OBSALT> Altitude of observatory (meters)
#<RA2000> Right ascension of object for epoch 2000.0 (hours)
#<DE2000> Declination of object for epoch 2000.0 (degrees)
#<JD> Julian date for the middle of exposure
#[DEBUG=] set keyword to get additional results for debugging
#
#OUTPUT:
#<CORRECTION> baricentric correction - correction for rotation of earth,
#   rotation of earth center about the eart-moon barycenter, eart-moon
#   barycenter about the center of the Sun.
#<HJD> Heliocentric Julian date for middle of exposure
#
#Algorithms used are taken from the IRAF task noao.astutils.rvcorrect
#and some procedures of the IDL Astrolib are used as well.
#Accuracy is about 0.5 seconds in time and about 1 m/s in velocity.
#
#History:
#written by Peter Mittermayer, Nov 8,2003
#2005-January-13   Kudryavtsev   Made more accurate calculation of the sideral time.
#                                Conformity with MIDAS compute/barycorr is checked.
#2005-June-20      Kochukhov Included precession of RA2000 and DEC2000 to current epoch

   
   #covert JD to Gregorian calendar date
   xjd = array(2400000.).astype(float) + jd
   year,month,day,ut=daycnv(xjd)
   
   #current epoch
   epoch = year + month / 12. + day / 365.
   
   #precess ra2000 and dec2000 to current epoch
   ra,dec=precess(ra2000*15., dec2000, 2000.0, epoch)
   #calculate heliocentric julian date
   hjd = array(helio_jd(jd, ra, dec)).astype(float)
   
   #DIURNAL VELOCITY (see IRAF task noao.astutil.rvcorrect)
   #convert geodetic latitude into geocentric latitude to correct
   #for rotation of earth
   dlat = -(11. * 60. + 32.743) * sin(2 * obs_lat / _radeg) + 1.1633 * sin(4 * obs_lat / _radeg) - 0.0026 * sin(6 * obs_lat / _radeg)
   lat = obs_lat + dlat / 3600
   
   #calculate distance of observer from earth center
   r = 6378160.0 * (0.998327073 + 0.001676438 * cos(2 * lat / _radeg) - 0.00000351 * cos(4 * lat / _radeg) + 0.000000008 * cos(6 * lat / _radeg)) + obs_alt
   
   #calculate rotational velocity (perpendicular to the radius vector) in km/s
   #23.934469591229 is the siderial day in hours for 1986
   v = 2. * pi * (r / 1000.) / (23.934469591229 * 3600.)
   
   #calculating local mean siderial time (see astronomical almanach)
   tu = (jd - 51545.0) / 36525
   gmst = 6.697374558 + ut + (236.555367908 * (jd - 51545.0) + 0.093104 * tu ** 2 - 6.2e-6 * tu ** 3) / 3600
   lmst = (gmst - obs_long / 15) % 24
   
   #projection of rotational velocity along the line of sight
   vdiurnal = v * cos(lat / _radeg) * cos(dec / _radeg) * sin((ra - lmst * 15) / _radeg)
   
   #BARICENTRIC and HELIOCENTRIC VELOCITIES
   vh,vb=baryvel(xjd, 0)
   
   #project to line of sight
   vbar = vb[0] * cos(dec / _radeg) * cos(ra / _radeg) + vb[1] * cos(dec / _radeg) * sin(ra / _radeg) + vb[2] * sin(dec / _radeg)
   vhel = vh[0] * cos(dec / _radeg) * cos(ra / _radeg) + vh[1] * cos(dec / _radeg) * sin(ra / _radeg) + vh[2] * sin(dec / _radeg)
   
   corr = (vdiurnal + vbar) #using baricentric velocity for correction
   
   if debug:   
      print ''
      print '----- HELCORR.PRO - DEBUG INFO - START ----'
      print '(obs_long,obs_lat,obs_alt) Observatory coordinates [deg,m]: ', obs_long, obs_lat, obs_alt
      print '(ra,dec) Object coordinates (for epoch 2000.0) [deg]: ', ra, dec
      print '(ut) Universal time (middle of exposure) [hrs]: ', ut#, format='(A,F20.12)'
      print '(jd) Julian date (middle of exposure) (JD-2400000): ', jd#, format='(A,F20.12)'
      print '(hjd) Heliocentric Julian date (middle of exposure) (HJD-2400000): ', hjd#, format='(A,F20.12)'
      print '(gmst) Greenwich mean siderial time [hrs]: ', gmst % 24
      print '(lmst) Local mean siderial time [hrs]: ', lmst
      print '(dlat) Latitude correction [deg]: ', dlat
      print '(lat) Geocentric latitude of observer [deg]: ', lat
      print '(r) Distance of observer from center of earth [m]: ', r
      print '(v) Rotational velocity of earth at the position of the observer [km/s]: ', v
      print '(vdiurnal) Projected earth rotation and earth-moon revolution [km/s]: ', vdiurnal
      print '(vbar) Baricentric velocity [km/s]: ', vbar
      print '(vhel) Heliocentric velocity [km/s]: ', vhel
      print '(corr) Vdiurnal+vbar [km/s]: ', corr#, format='(A,F12.9)'
      print '----- HELCORR.PRO - DEBUG INFO - END -----'
      print ''
   
   
   return (corr, hjd)


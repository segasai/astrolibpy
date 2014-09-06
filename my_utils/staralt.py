# -*- coding: utf-8 -*-
# Copyright (C) 2014 Sergey Koposov koposov@ast.cam.ac.uk
#
# This file is part of astrolibpy
#
#    astrolibpy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#   astrolibpy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with astrolibpy.  If not, see <http://www.gnu.org/licenses/>.

import ephem,numpy as np

import sphdist

def getalt(ra, dec, yr, mon, day, hr, minu, lon='-111:35:59', lat='31:57:12',
		elev=2000, retSun=False, retDistMoon=False):
	#computes the altitude in degrees of a given object at the given utc time
	#ra dec in degress	
	#yr mon day hr minu in utc
	# lon lat are in degrees and lon is positive to the East
	
	obs = ephem.Observer();
	obs.lon = lon # longi tude 
	obs.lat = lat #latitude
	obs.elevation=elev;
	fb = ephem.FixedBody();
	fbSun = ephem.Sun();
	fbMoon = ephem.Moon();
	fb._ra=np.deg2rad(ra);
	fb._dec=np.deg2rad(dec)
	obs.date=ephem.Date('%d/%02d/%d %d:%d:0'%(yr,mon,day,hr,minu))
	fb.compute(obs);
	fbSun.compute(obs);
	fbMoon.compute(obs);
	alt = np.rad2deg(1*fb.alt)
	az = np.rad2deg(1*fb.az);
	altSun = np.rad2deg(1*(fbSun.alt))
	altMoon = np.rad2deg(1*fbMoon.alt)
	azMoon = np.rad2deg(1*fbMoon.az);
	distMoon = sphdist.sphdist(az,alt,azMoon,altMoon)
	ret=[alt]
	if retSun:
		ret.append(altSun)
	if retDistMoon:
		ret.append(distMoon)
	return tuple(ret)
	

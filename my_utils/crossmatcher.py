import sqlutilpy
import numpy as np


def doit(tabname, ra, dec, colstring, radeccols=('ra', 'dec'), rad=1.,
		 extra=None, yourradeccols=('ra', 'dec'), host=None,
		 db=None,user=None,password=None):
	"""
	Performs the nearest neighbor crossmatch within specified radius in arcseconds
	with the remote table the input is given by ra,dec columns 
	
	colstring option is for the comma sepated list of columns that you want to
	retrieve	
	radeccols option allows you to specify the name of ra,dec columns 
	on the remote table. 
	rad options is the cross-match radius in arcseconds
	extra allows you to specify and additional SQL Where condition 
	Example:
	ra = np.arange(10)
	dec = np.arange(10)+5
	gmag,rmag= crossmatcher.doit('sdssdr9.phototag', ra,dec,
			'psfmag_g,psfmag_r', rad=2.)
	
	"""
	racol, deccol = radeccols
	if extra is None:
		extra = 'true'
	your_ra, your_dec = yourradeccols  
	RES = sqlutilpy.local_join(str.format("""
		select {colstring} from mytable as m 
			left join lateral (select * from {tabname} as s 
					where 
					q3c_join(m.{your_ra}, m.{your_dec}, 
							s.{racol}, s.{deccol}, {rad}/3600.) 
					and {extra}
					order by 
			q3c_dist(m.{your_ra}, m.{your_dec}, 
						s.{racol},s.{deccol}) asc limit 1) as tt 
				on true order by xid """, **locals()), 'mytable',
		(ra, dec, np.arange(len(ra))), (your_ra, your_dec, 'xid'),
		preamb='set enable_seqscan to off; set enable_mergejoin to off; set enable_hashjoin to off;',
		host=host,db=db,user=user,password=password)
	return RES

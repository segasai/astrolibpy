import sqlutilpy
import numpy as np


def doit(tabname,
         ra,
         dec,
         colstring,
         radeccols=('ra', 'dec'),
         rad=1.,
         extra=None,
         yourradeccols=('ra', 'dec'),
         tab_alias='tt',
         host=None,
         port=None,
         db=None,
         user=None,
         password=None,
         asDict=False,
         conn=None,
         preamb=None):
    """
    Performs the nearest neighbor crossmatch within specified radius
    with the remote table the input is given by ra,dec columns

    Parameters:
    -----------

    tabname: string
        The name of the table to crossmatch against
    ra: numpy
        The numpy array with right ascension
    dec : numpy
        The numpy array with declinations
    colstring: string
        The comma sepated list of columns that you want to retrieve
    radeccols: tuple (optional)
        The tuple of two strings with the name of ra,dec columns in
        the remote table. Default ('ra','dec')
    rad: float (optional)
        The cross-match radius in arcseconds (default 1)
    extra: string allows you to specify and additional SQL Where condition
    tab_alias: string
        The alias for the table that you are crossmatching against, so you can
        refer to its columns
    host,port,db,user,password: string
        The connection parameters to the DB if needed
    asDict: bool
        if True instead of returning a tuple a dictionary is returned
    preamb: string (optional)
        additional commands to be executed before the query
    conn: connection object(optional)
        An explicit connection to the DB can be provided

    Example:
    --------
    > ra = np.arange(10)
    > dec = np.arange(10)+5
    > gmag,rmag= crossmatcher.doit('sdssdr9.phototag', ra,dec,
           'psfmag_g,psfmag_r', rad=2.)

    """
    racol, deccol = radeccols
    if extra is None:
        extra = 'true'
    your_ra, your_dec = yourradeccols
    preamb = '' or preamb
    RES = sqlutilpy.local_join(
        str.format(
            """
            select {colstring} from
                  ( select * from mytable order by
            q3c_ang2ipix({your_ra},{your_dec})
                  ) as m
            left join lateral (select * from {tabname} as s where
            q3c_join(m.{your_ra}, m.{your_dec},
            s.{racol}, s.{deccol}, {rad}/3600.) and {extra}
                order by q3c_dist(m.{your_ra}, m.{your_dec},
                s.{racol},s.{deccol}) asc limit 1) as {tab_alias}
            on true order by xid """, **locals()),
        'mytable', (ra, dec, np.arange(len(ra))), (your_ra, your_dec, 'xid'),
        preamb=('set enable_seqscan to off;' + 'set enable_mergejoin to off;' +
                'set enable_hashjoin to off;' + (preamb or '')),
        host=host,
        db=db,
        user=user,
        port=(port or 5432),
        password=password,
        asDict=asDict,
        conn=conn)
    return RES

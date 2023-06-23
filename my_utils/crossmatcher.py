import sqlutilpy
import numpy as np


def doit(tabname,
         ra,
         dec,
         colstring,
         radeccols=('ra', 'dec'),
         pmradeccols=('pmra', 'pmdec'),
         epoch=None,
         rad=1.,
         extra=None,
         epochcol=None,
         max_pm_offset=30,
         yourradeccols=('ra', 'dec'),
         yourepochcol='epoch',
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
    epoch: numpy (optional)
        If specified we try to use proper motions for the xmatch
    radeccols: tuple (optional)
        The tuple of two strings with the name of ra,dec columns in
        the remote table. Default ('ra','dec')
    pmradeccols: tuple (optional)
        The tuple of proper motion columns in the table
    epochcol: str
        The name of the column with epoch
    rad: float (optional)
        The cross-match radius in arcseconds (default 1)
    max_pm_offset: float (optional)
        The maximum offset in arcsec allowed. If proper motions
        are used
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
    kw = dict(
        preamb=('set enable_seqscan to off;' + 'set enable_mergejoin to off;' +
                'set enable_hashjoin to off;' + (preamb or '')),
        host=host,
        db=db,
        user=user,
        port=(port or 5432),
        password=password,
        asDict=asDict,
        conn=conn)

    if epoch is None:
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
            on true order by xid """, **locals()), 'mytable',
            (ra, dec, np.arange(len(ra))), (your_ra, your_dec, 'xid'), **kw)
    else:
        try:
            nep = len(epoch)
        except TypeError:
            epoch = np.zeros(len(ra)) + epoch
        else:
            if nep != len(ra):
                raise Exception(
                    'length of the epoch must be equal to length of positions')
        maxap = max_pm_offset + rad
        pmracol, pmdeccol = pmradeccols
        dist_str = f'''q3c_dist_pm(
                s.{racol},s.{deccol}, s.{pmracol}, s.{pmdeccol}, 1, 
                s.{epochcol}, m.{your_ra}, m.{your_dec}, 
                m.{yourepochcol})'''
        RES = sqlutilpy.local_join(
            str.format(
                """
            select {colstring} from
                  ( select * from mytable order by
            q3c_ang2ipix({your_ra},{your_dec})
                  ) as m
            left join lateral (select * from {tabname} as s where
            q3c_join(m.{your_ra}, m.{your_dec},
            s.{racol}, s.{deccol}, {maxap}/3600.) and {extra}
                and {dist_str}<{rad}/3600.
                order by {dist_str} asc limit 1) as {tab_alias}
            on true order by xid """,
                **locals()), 'mytable', (ra, dec, np.arange(len(ra)), epoch),
            (your_ra, your_dec, 'xid', yourepochcol), **kw)

    return RES


def doit_by_key(tabname,
                keys,
                colstring,
                key_col=None,
                rad=1.,
                extra=None,
                yourkeycol='id',
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
    Performs the crossmatch by id
    with the remote table 

    Parameters:
    -----------

    tabname: string
        The name of the table to crossmatch against
    keys: numpy
        The numpy array with ids that will be used for matching
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
    if extra is None:
        extra = 'true'
    preamb = '' or preamb
    RES = sqlutilpy.local_join(
        str.format(
            """
            select {colstring} from
                  ( select * from mytable
                  ) as m
            left join lateral (select * from {tabname} as s where
            m.id=s.{key_col} and {extra} limit 1) as {tab_alias}
            on true order by xid """, **locals()),
        'mytable', (keys, np.arange(len(keys))), ('id', 'xid'),
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

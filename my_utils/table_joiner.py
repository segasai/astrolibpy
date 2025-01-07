import duckdb
import numpy as np
import pandas as pd
import astropy.table as atpy


def order_convert(x):
    if x.dtype.byteorder not in ('=', '|'):
        x = x.byteswap()
        return x.view(x.data.dtype.newbyteorder())
    return x


def check_type(arr):
    # Check if the array is of string type
    if np.issubdtype(arr.dtype, np.str_):
        return 'string'
    if np.issubdtype(arr.dtype, np.bytes_):
        return 'string'
    elif np.issubdtype(arr.dtype, np.integer):
        return 'int'
    if np.issubdtype(arr.dtype, np.floating):
        return 'float'
    raise RuntimeError('unsupported type')


def join_to_left(T1,
                 T2,
                 key=None,
                 fill_value=None,
                 key_left=None,
                 key_right=None,
                 fill_values={
                     'string': '',
                     'int': -9999,
                     'float': np.nan
                 }):
    """
    Return Table 2 matched to T1 using key.
    
    In the case of multiple matches a random one is returned

    """
    rowid1 = np.arange(len(T1))
    assert ((key is not None)
            or (key_left is not None and key_right is not None))
    if key is not None:
        key_left, key_right = key, key

    df1 = pd.DataFrame({'rowid': rowid1, 'key': order_convert(T1[key_left])})
    rowid2 = np.arange(len(T2))
    df2 = pd.DataFrame({'rowid': rowid2, 'key': order_convert(T2[key_right])})
    duckdb.register("tab1", df1)
    duckdb.register("tab2", df2)
    R = duckdb.sql('''
    select distinct on (tab1.rowid) tab1.rowid as r1, tab2.rowid  as r2
    from tab1, tab2 where
    tab1.key = tab2.key order by tab1.rowid''').fetchnumpy()
    xrowid1, xrowid2 = R['r1'], R['r2']
    ret = {}
    n1 = len(T1)
    for k in T2.columns:
        cura = T2[k][xrowid2]
        cur_fill = fill_values[check_type(cura)]
        cur_ret = np.zeros(n1, dtype=cura.dtype)
        cur_ret[:] = cur_fill
        cur_ret[xrowid1] = cura
        ret[k] = cur_ret
    return atpy.Table(ret)

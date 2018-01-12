import astropy.table as atpy2
import pandas
import sqlalchemy
import sqlutil

def printSchema(fname, host, db, user='koposov', tabname='newtab'):
    # print the sql schema for a given table 
    tab = atpy2.Table().read(fname)
    df = tab.filled().to_pandas()
    conn = sqlutil.getConnection(driver='psycopg2', host=host, db=db, user=user)
    AL = sqlalchemy.create_engine('postgresql://', creator=(lambda : conn))
    print (pandas.io.sql.get_schema(df,tabname, con=AL))


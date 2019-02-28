"""
Accumulate all experiment data into a set of linked tables.
"""
import os, io, time
from datetime import datetime
import numpy as np

import sqlalchemy
from distutils.version import LooseVersion
if LooseVersion(sqlalchemy.__version__) < '1.2':
    raise Exception('requires at least sqlalchemy 1.2')

from sqlalchemy import create_engine, Column, Integer, String, Boolean, Float, Date, DateTime, LargeBinary, ForeignKey, or_, and_
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship, deferred, sessionmaker, aliased
from sqlalchemy.types import TypeDecorator
from sqlalchemy.sql.expression import func

from .. import config

# database version should be incremented whenever the schema has changed
db_version = 11
db_name = '{database}_{version}'.format(database=config.synphys_db, version=db_version)
db_address_ro = '{host}/{database}'.format(host=config.synphys_db_host, database=db_name)
if config.synphys_db_host_rw is None:
    db_address_rw = None
else:
    db_address_rw = '{host}/{database}'.format(host=config.synphys_db_host_rw, database=db_name)

default_sample_rate = 20000


_sample_rate_str = '%dkHz' % (default_sample_rate // 1000)


class TableGroup(object):
    """Class used to manage a group of tables that act as a single unit--tables in a group
    are always created and deleted together.
    """
    schemas = {}
    
    def __init__(self):
        self.mappings = {}
        self.create_mappings()

    def __getitem__(self, item):
        return self.mappings[item]

    def create_mappings(self):
        for k,schema in self.schemas.items():
            self.mappings[k] = generate_mapping(k, schema)

    def drop_tables(self):
        global engine_rw
        for k in self.schemas:
            if k in engine_rw.table_names():
                self[k].__table__.drop(bind=engine_rw)

    def create_tables(self):
        global engine_rw, engine_ro
        if engine_rw is None:
            for k in self.schemas:
                if k not in engine_ro.table_names():
                    raise Exception("Table %s not found in database %s" % (k, db_address_ro))
            return
        for k in self.schemas:
            if k not in engine_rw.table_names():
                self[k].__table__.create(bind=engine_rw)





#----------- define ORM classes -------------

ORMBase = declarative_base()

class NDArray(TypeDecorator):
    """For marshalling arrays in/out of binary DB fields.
    """
    impl = LargeBinary
    
    def process_bind_param(self, value, dialect):
        if value is None:
            return b'' 
        buf = io.BytesIO()
        np.save(buf, value, allow_pickle=False)
        return buf.getvalue()
        
    def process_result_value(self, value, dialect):
        if value == b'':
            return None
        buf = io.BytesIO(value)
        return np.load(buf, allow_pickle=False)


class FloatType(TypeDecorator):
    """For marshalling float types (including numpy).
    """
    impl = Float
    
    def process_bind_param(self, value, dialect):
        if value is None:
            return None
        return float(value)
        
    #def process_result_value(self, value, dialect):
        #buf = io.BytesIO(value)
        #return np.load(buf, allow_pickle=False)


_coltypes = {
    'int': Integer,
    'float': FloatType,
    'bool': Boolean,
    'str': String,
    'date': Date,
    'datetime': DateTime,
    'array': NDArray,
    'object': JSONB,
}


def generate_mapping(table, schema, base=None):
    """Generate an ORM mapping class from an entry in table_schemas.
    """
    name = table.capitalize()
    table_args = {}
    if isinstance(schema[0], str):
        table_args['comment'] = schema[0]
        schema = schema[1:]
    
    props = {
        '__tablename__': table,
        '__table_args__': table_args,
        'id': Column(Integer, primary_key=True),
    }
    for column in schema:
        colname, coltype = column[:2]
        kwds = {} if len(column) < 4 else column[3]
        kwds['comment'] = None if len(column) < 3 else column[2]
        defer_col = kwds.pop('deferred', False)
        ondelete = kwds.pop('ondelete', None)

        if coltype not in _coltypes:
            if not coltype.endswith('.id'):
                raise ValueError("Unrecognized column type %s" % coltype)
            props[colname] = Column(Integer, ForeignKey(coltype, ondelete=ondelete), **kwds)
        else:
            ctyp = _coltypes[coltype]
            props[colname] = Column(ctyp, **kwds)

        if defer_col:
            props[colname] = deferred(props[colname])

    props['time_created'] = Column(DateTime, default=func.now())
    props['time_modified'] = Column(DateTime, onupdate=func.current_timestamp())
    props['meta'] = Column(JSONB)

    if base is None:
        return type(name, (ORMBase,), props)
    else:
        def init(self, *args, **kwds):
            base.__init__(self)
            ORMBase.__init__(self, *args, **kwds)
        props['__init__'] = init  # doesn't work?
        return type(name, (base,ORMBase), props)


def _generate_mapping(table, base=None):
    return generate_mapping(table, table_schemas[table], base=base)




#-------------- initial DB access ----------------
engine_ro = None
engine_rw = None
engine_pid = None  # pid of process that created this engine. 
def init_engine():
    global engine_ro, engine_rw, engine_pid
    dispose_engines()
    
    engine_ro = create_engine(db_address_ro, pool_size=10, max_overflow=40, isolation_level='AUTOCOMMIT')
    if db_address_rw is not None:
        engine_rw = create_engine(db_address_rw, pool_size=10, max_overflow=40)
    engine_pid = os.getpid()


def dispose_engines():
    global engine_ro, engine_rw, engine_pid
    if engine_ro is not None:
        engine_ro.dispose()
        engine_ro = None
    if engine_rw is not None:
        engine_rw.dispose()
        engine_rw = None
    engine_pid = None    

init_engine()


_sessionmaker_ro = None
_sessionmaker_rw = None
# external users should create sessions from here.
def Session(readonly=True):
    """Create and return a new database Session instance.
    
    If readonly is True, then the session is created using read-only credentials and has autocommit enabled.
    This prevents idle-in-transaction timeouts that occur when GUI analysis tools would otherwise leave transactions
    open after each request.
    """
    global _sessionmaker_ro, _sessionmaker_rw, engine_ro, engine_rw, engine_pid
    if os.getpid() != engine_pid:
        # In forked processes, we need to re-initialize the engine before
        # creating a new session, otherwise child processes will
        # inherit and muck with the same connections. See:
        # http://docs.sqlalchemy.org/en/rel_1_0/faq/connections.html#how-do-i-use-engines-connections-sessions-with-python-multiprocessing-or-os-fork
        # if engine_pid is not None:
        #     print("Making new session for subprocess %d != %d" % (os.getpid(), engine_pid))
        init_engine()
        _sessionmaker_ro = None
    
    if _sessionmaker_ro is None:
        _sessionmaker_ro = sessionmaker(bind=engine_ro)
    if _sessionmaker_rw is None and engine_rw is not None:    
        _sessionmaker_rw = sessionmaker(bind=engine_rw)
        
    if readonly:
        return _sessionmaker_ro()
    else:
        return _sessionmaker_rw()


create_all_mappings()



def reset_db():
    """Drop the existing synphys database and initialize a new one.
    """
    global engine_rw
    
    pg_engine = create_engine(config.synphys_db_host_rw + '/postgres')
    with pg_engine.begin() as conn:
        conn.connection.set_isolation_level(0)
        try:
            conn.execute('drop database %s' % db_name)
        except sqlalchemy.exc.ProgrammingError as err:
            if 'does not exist' not in err.message:
                raise

        conn.execute('create database %s' % db_name)

    # reconnect to DB
    init_engine()

    # Grant readonly permissions
    ro_user = config.synphys_db_readonly_user
    if ro_user is not None:
        with engine_rw.begin() as conn:
            conn.execute('ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO %s;' % ro_user)
            # should only be needed if there are already tables present
            #conn.execute('GRANT SELECT ON ALL TABLES IN SCHEMA public TO %s;' % ro_user)
    
    # Create all tables
    global ORMBase
    ORMBase = declarative_base()
    create_all_mappings()
    ORMBase.metadata.create_all(engine_rw)


def vacuum(tables=None):
    """Cleans up database and analyzes table statistics in order to improve query planning.
    Should be run after any significant changes to the database.
    """
    global engine_rw
    with engine_rw.begin() as conn:
        conn.connection.set_isolation_level(0)
        if tables is None:
            conn.execute('vacuum analyze')
        else:
            for table in tables:
                conn.execute('vacuum analyze %s' % table)


def default_session(fn):
    def wrap_with_session(*args, **kwds):
        close = False
        if kwds.get('session', None) is None:
            kwds['session'] = Session(readonly=True)
            close = True
        try:
            ret = fn(*args, **kwds)
            return ret
        finally:
            if close:
                kwds['session'].close()
    return wrap_with_session    


@default_session
def slice_from_timestamp(ts, session=None):
    slices = session.query(Slice).filter(Slice.acq_timestamp==ts).all()
    if len(slices) == 0:
        raise KeyError("No slice found for timestamp %0.3f" % ts)
    elif len(slices) > 1:
        raise KeyError("Multiple slices found for timestamp %0.3f" % ts)
    
    return slices[0]


@default_session
def experiment_from_timestamp(ts, session=None):
    expts = session.query(Experiment).filter(Experiment.acq_timestamp==ts).all()
    if len(expts) == 0:
        # For backward compatibility, check for timestamp truncated to 2 decimal places
        for expt in session.query(Experiment).all():
            if abs((expt.acq_timestamp - ts)) < 0.01:
                return expt
        
        raise KeyError("No experiment found for timestamp %0.3f" % ts)
    elif len(expts) > 1:
        raise RuntimeError("Multiple experiments found for timestamp %0.3f" % ts)
    
    return expts[0]


@default_session
def list_experiments(session=None):
    return session.query(Experiment).all()


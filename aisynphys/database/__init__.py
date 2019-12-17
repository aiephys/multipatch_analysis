from .. import config
from .database import Database
from .synphys_database import SynphysDatabase


class NoDatabase:
    """Raises an exception on first attempt to access an attribute.
    """
    def __init__(self, msg):
        self.msg = msg
    def __getattr__(self, attr):
        raise Exception(self.msg)


# initialize a default database connection if configured or requested via CLI
if config.synphys_db_host is None:
    default_db = NoDatabase("No database was specified in config.synphys_db_host or with CLI flags --db-version or --db-host")
else:
    if config.synphys_db_host.startswith('postgres'):
        default_db_name = '{database}_{version}'.format(database=config.synphys_db, version=SynphysDatabase.schema_version)
    else:
        default_db_name = config.synphys_db
    default_db = SynphysDatabase(config.synphys_db_host, config.synphys_db_host_rw, default_db_name)


def dispose_all_engines():
    """Dispose all engines across all Database instances.
    
    This function should be called before forking.
    """
    Database.dispose_all_engines()

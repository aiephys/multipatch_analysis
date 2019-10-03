from .. import config
from .database import Database
from .synphys_database import SynphysDatabase

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

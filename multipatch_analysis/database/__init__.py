from .. import config
from .database import Database, db_version, aliased, or_, and_, default_sample_rate
from .synphys_database import SynphysDatabase

    
default_db_name = '{database}_{version}'.format(database=config.synphys_db, version=db_version)
default_db = SynphysDatabase(config.synphys_db_host, config.synphys_db_host_rw, default_db_name)


def dispose_all_engines():
    """Dispose all engines across all Database instances.
    
    This function should be called before forking.
    """
    Database.dispose_all_engines()

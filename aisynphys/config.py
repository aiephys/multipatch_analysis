"""
Site-specific configuration parameters.

Local variables in this module are overwritten by the contents of config.yml

"""

import os, sys, yaml, argparse

# default cache path in user's home dir
cache_path = os.path.join(os.path.expanduser('~'), 'ai_synphys_cache')

# Parameters for the DB connection provided by aisynphys.database.default_db
# For sqlite files:
#    synphys_db_host = "sqlite:///"
#    synphys_db = "path/to/database.sqlite"
# For postgres
#    synphys_db_host = "postgresql://user:password@hostname"
#    synphys_db = "database_name"
synphys_db_host = None
synphys_db = None


# utility config, not meant for external use
synphys_data = None  # location of data repo network storage
synphys_db_host_rw = None  # rw access to postgres / sqlite DB
synphys_db_readonly_user = "readonly"  # readonly postgres username assigned whrn creating db/tables
lims_address = None
rig_name = None
n_headstages = 8
rig_data_paths = {}
known_addrs = {}


configfile = os.path.join(os.path.dirname(__file__), '..', 'config.yml')

if os.path.isfile(configfile):
    if hasattr(yaml, 'FullLoader'):
        # pyyaml new API
        config = yaml.load(open(configfile, 'rb'), Loader=yaml.FullLoader)
    else:
        # pyyaml old API
        config = yaml.load(open(configfile, 'rb'))

    if config is None:
        config = {}

    for k,v in config.items():
        locals()[k] = v


# intercept specific command line args
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--db-version', default=None, dest='db_version')
parser.add_argument('--db-host', default=None, dest='db_host')
parser.add_argument('--database', default=None)

args, unknown_args = parser.parse_known_args()
sys.argv = sys.argv[:1] + unknown_args

if args.db_version is not None:
    from .synphys_cache import get_db_path
    sqlite_file = get_db_path(args.db_version)
    synphys_db_host = "sqlite:///"
    synphys_db_host_rw = None
    synphys_db = sqlite_file
    
else:
    if args.db_host is not None:
        synphys_db_host = args.db_host
    if args.database is not None:
        synphys_db = args.database
        

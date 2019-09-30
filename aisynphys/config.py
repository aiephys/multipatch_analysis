"""
Site-specific configuration parameters.

Local variables in this module are overwritten by the contents of config.yml

"""

import os, sys, yaml

# default cache path in user's home dir
cache_path = os.path.join(os.path.expanduser('~'), 'ai_synphys_cache')

synphys_db_host = "sqlite:///"
synphys_db_host_rw = None
synphys_db = "synphys"
synphys_db_readonly_user = "readonly"
synphys_data = None
lims_address = None
grow_cache = False
rig_name = None
n_headstages = 8
raw_data_paths = []
rig_data_paths = {}
known_addrs = {}
import_old_data_on_submission = False


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
ignored_args = [sys.argv[0]]
for arg in sys.argv[1:]:
    if arg.startswith('--database='):
        synphys_db = arg[11:]
    elif arg.startswith('--db-host='):
        synphys_db_host = arg[10:]
    else:
        ignored_args.append(arg)
sys.argv = ignored_args        

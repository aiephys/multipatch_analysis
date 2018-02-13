"""
Site-specific configuration parameters.

Local variables in this module are overwritten by the contents of config.yml

"""

import os, yaml


synphys_db_host = None
synphys_db = "synphys"
synphys_db_readonly_user = None
synphys_data = None
cache_path = "cache"
rig_name = None
n_headstages = 8
raw_data_paths = []
summary_files = []


template = """
synphys_db_host: "postgresql://postgres:xxxxx@10.128.38.98"
synphys_db: "synphys"
synphys_data: "/path/to/server/synphys_data"
cache_path: "cache"
rig_name: 'MP_'
n_headstages: 8
raw_data_paths:
    - '/raw/data/path/1'
    - '/raw/data/path/2'
summary_files:
    - '/path/to/old/connectivity_summary'
    - '/path/to/old/connectivity_summary'    
"""

configfile = os.path.join(os.path.dirname(__file__), '..', 'config.yml')
if not os.path.isfile(configfile):
    open(configfile, 'wb').write(template)

config = yaml.load(open(configfile, 'rb'))

for k,v in config.items():
    locals()[k] = v




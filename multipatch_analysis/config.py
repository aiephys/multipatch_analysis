import os, yaml

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
"""

configfile = os.path.join(os.path.dirname(__file__), 'config.yml')
if not os.path.isfile(configfile):
    open(configfile, 'wb').write(template)

config = yaml.load(open(configfile, 'rb'))

for k,v in config.items():
    locals()[k] = v




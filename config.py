import os, yaml

template = """
synphys_db: "postgresql://user:password@server/synphys"
synphys_data: "/path/to/data"
cache_path: "cache"
rig_name: 'MP2'
n_headstages: 8
"""

configfile = 'config.yml'
if not os.path.isfile(configfile):
    open(configfile, 'wb').write(template)

config = yaml.load(open('config.yml', 'rb'))

for k,v in config.items():
    locals()[k] = v




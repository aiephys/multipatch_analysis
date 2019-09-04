"""
Site-specific configuration parameters.

Local variables in this module are overwritten by the contents of config.yml

"""

import os, sys, yaml


synphys_db_host = None
synphys_db_host_rw = None
synphys_db = "synphys"
synphys_db_readonly_user = "readonly"
synphys_data = None
lims_address = None
cache_path = "cache"
grow_cache = False
rig_name = None
n_headstages = 8
raw_data_paths = []
rig_data_paths = {}
known_addrs = {}
import_old_data_on_submission = False

template = r"""
synphys_db_host: "sqlite:///"
synphys_db: "synphys"

cache_path: "E:\\multipatch_analysis_cache"
grow_cache: true

editor_command: '"C:\\Program Files\\Sublime Text 2\\sublime_text.exe" "{file}"'
browser_command: '"C:\\Program Files (x86)\\Mozilla Firefox\\firefox.exe" {url}'
        
"""

configfile = os.path.join(os.path.dirname(__file__), '..', 'config.yml')
if not os.path.isfile(configfile):
    open(configfile, 'wb').write(template.encode('utf8'))

if hasattr(yaml, 'FullLoader'):
    # pyyaml new API
    config = yaml.load(open(configfile, 'rb'), Loader=yaml.FullLoader)
else:
    # pyyaml old API
    config = yaml.load(open(configfile, 'rb'))

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



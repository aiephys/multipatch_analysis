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


template = r"""
# synphys database
synphys_db_host: "postgresql://readonly:readonly@10.128.36.109"
synphys_db: "synphys"

# path to synphys network storage
synphys_data: "N:\\"

cache_path: "E:\\multipatch_analysis_cache"
rig_name: 'MP_'
n_headstages: 8

editor_command: '"C:\\Program Files\\Sublime Text 2\\sublime_text.exe" "{file}"'
browser_command: '"C:\\Program Files (x86)\\Mozilla Firefox\\firefox.exe" {url}'

# local paths to data sources
rig_data_paths:
    mp1:
        - primary: /path/to/mp1/primary_data_1
          archive: /path/to/mp1/data_1_archive
          backup:  /path/to/mp1/data_1_backup
        - primary: /path/to/mp1/primary_data_2
          archive: /path/to/mp1/data_2_archive
          backup:  /path/to/mp1/data_2_backup


# directories to be synchronized nightly
backup_paths:
    rig_data:
        source: "D:\\"
        dest: "E:\\archive"
        archive_deleted: false
    system_drive:
        source: "C:\\"
        dest: "E:\\C_backup"
        archive_deleted: true
        
"""

configfile = os.path.join(os.path.dirname(__file__), '..', 'config.yml')
if not os.path.isfile(configfile):
    open(configfile, 'wb').write(template)

config = yaml.load(open(configfile, 'rb'))

for k,v in config.items():
    locals()[k] = v




import os, yaml

template = """
synphys_db_host: "postgresql://postgres:xxxxx@10.128.38.98"
synphys_db: "synphys"
synphys_data: "/path/to/server/synphys_data"
cache_path: "cache"
rig_name: 'MP_'
n_headstages: 8
raw_data_paths:
    - '/home/luke/mnt/backup_server/MP1_backup/D_drive/Steph/'
    - '/home/luke/mnt/backup_server/MP2_backup/D_drive/data/Pasha/V1'
    - '/home/luke/mnt/backup_server/MP2_backup/D_drive/data/Pasha/Human'
    - '/home/luke/mnt/backup_server/MP3_backup/D_drive/data/Alex/V1/'
    - '/home/luke/mnt/backup_server/MP3_backup/D_drive/data/Alex/Human/'
    - '/home/luke/mnt/backup_server/MP3_backup/version_backups/data/Alex/V1/'
    - '/home/luke/mnt/backup_server/MP3_backup/version_backups/data/Alex/Human/'

"""

configfile = os.path.join(os.path.dirname(__file__), 'config.yml')
if not os.path.isfile(configfile):
    open(configfile, 'wb').write(template)

config = yaml.load(open(configfile, 'rb'))

for k,v in config.items():
    locals()[k] = v




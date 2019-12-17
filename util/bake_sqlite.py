import os, sys, datetime
from aisynphys.database import default_db as db

date = datetime.datetime.today().strftime("%Y-%m-%d")

db_files = {
    'small': "db_bakes/synphys_%s_small.sqlite" % date,
    'medium': "db_bakes/synphys_current_medium.sqlite",
    'full': "db_bakes/synphys_current_full.sqlite",
}

skip_tables = {}
skip_tables['full'] = []
skip_tables['medium'] = skip_tables['full'] + []
skip_tables['small'] = skip_tables['medium'] + [
    'synapse_prediction',
    'pulse_response_strength',
    'baseline_response_strength',
    'pulse_response_fit',
    'pipeline',
    'pulse_response',
    'baseline',
    'stim_pulse',
    'stim_spike',
    'recording',
    'patch_clamp_recording',
    'multi_patch_probe',
    'test_pulse',
    'sync_rec',    
]

skip_columns = {}
skip_columns['full'] = {}
skip_columns['medium'] = skip_columns['full'].copy()
skip_columns['medium'].update({
    'pulse_response': ['data'],
    'baseline': ['data'],
    'stim_pulse': ['data'],
})
skip_columns['small'] = skip_columns['medium'].copy()


versions = sys.argv[1:]

if len(versions) == 0:
    versions = list(db_files.keys())


for version in versions:
    filename = db_files[version]
    print("========== Cloning %s DB %s =============" % (version, filename))
    if os.path.exists(filename):
        os.remove(filename)
    db.bake_sqlite(filename, skip_tables=skip_tables[version], skip_columns=skip_columns[version])

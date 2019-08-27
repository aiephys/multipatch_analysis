import os, datetime
from multipatch_analysis.database import default_db as db

medium_table_skips = [
]

medium_column_skips = {
    'pulse_response': ['data'],
    'baseline': ['data'],
    'stim_pulse': ['data'],
}

small_table_skips = medium_table_skips + [
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

small_column_skips = medium_column_skips.copy()
small_column_skips.update({
    
})

date = datetime.datetime.today().strftime("%Y-%m-%d")

small_file = "synphys_%s_small.sqlite" % date
print("========== Cloning small DB %s =============" % small_file)
if os.path.exists(small_file):
    os.remove(small_file)
db.bake_sqlite(small_file, skip_tables=small_table_skips, skip_columns=small_column_skips)

med_file = "synphys_current_medium.sqlite"
print("========== Cloning medium DB %s =============" % med_file)
if os.path.exists(med_file):
    os.remove(med_file)
db.bake_sqlite(med_file, skip_tables=medium_table_skips, skip_columns=medium_column_skips)

full_file = "synphys_current_full.sqlite"
print("========== Cloning full DB %s =============" % full_file)
if os.path.exists(full_file):
    os.remove(full_file)
db.bake_sqlite(full_file)



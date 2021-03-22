import sys, argparse
from aisynphys import config
from aisynphys.database import default_db as db

parser = argparse.ArgumentParser(parents=[config.parser])
parser.add_argument('experiment_id', type=str)
parser.add_argument('pre_cell_id', type=str, nargs='?')
parser.add_argument('post_cell_id', type=str, nargs='?')

args = parser.parse_args()

cell_ids = [cid for cid in (args.pre_cell_id, args.post_cell_id) if cid is not None]

expt = db.experiment_from_ext_id(args.experiment_id)
pipeline_jobs = {job.module_name:job for job in db.query(db.Pipeline).filter(db.Pipeline.job_id==args.experiment_id).all()}


def print_attrs(obj, attrs, indent=0, prefix=""):
    maxlen = max([len(a) for a in attrs])
    for attr in attrs:
        if isinstance(attr, tuple):
            attr, fmt, fn = attr
        else:
            fmt = str(maxlen) + "s"
            fn = str
        print(" "*indent + ("{prefix}{attr:" + fmt + "}  :  {val:s}").format(attr=attr, val=fn(getattr(obj, attr)), prefix=prefix))


def print_heading(txt, level=0):
    if level == 0:
        print("\n========================================")
        print(f"    {txt}")
        print("========================================")
    elif level == 1:
        print(f"\n    {txt}")
        print("    ----------------------------------------")
    elif level == 2:
        print(f"\n        --- {txt}")

print_heading(f"Experiment {expt.ext_id}")

print_attrs(expt, [
    'project_name', 
    'date', 
    'target_region', 
    'target_temperature',
    'internal', 
    'acsf', 
    'rig_name', 
    'operator_name', 
    'ephys_file', 
    'original_path', 
])


print_heading(f"Slice {expt.slice.ext_id}")

print_attrs(expt.slice, [
    'species',
    'age',
    'sex',
    'quality',
    'date_of_birth',
    'genotype',
    'hemisphere',
    'orientation',
    'surface',
])


print_heading(f"Cells")

for cell in expt.cell_list:
    if len(cell_ids) > 0 and cell.ext_id not in cell_ids:
        continue
    print_heading(f"Cell {cell.ext_id}", level=1)

    print_attrs(cell, [
        "cell_class",
        "cell_class_nonsynaptic",
        "cre_type",
        "depth",
        "target_layer",
        "meta",
    ], indent=4)

    if cell.cortical_location is None:
        print("    cortical_location: None")
    else:
        print_attrs(cell.cortical_location, ['cortical_layer'], indent=4, prefix="cortical_location.")

print_heading(f"Pairs")

pair_ids = sorted(list(expt.pairs.keys()))
for pair_id in pair_ids:
    if len(cell_ids) == 1 and (cell_ids[0] not in pair_id):
        continue
    if len(cell_ids) == 2 and tuple(cell_ids) != pair_id:
        continue
    pair = expt.pairs[pair_id]
    print_heading(f"Pair {expt.ext_id} {pair_id[0]} {pair_id[1]}", level=1)

    print_attrs(pair, [
        "has_synapse",
        "has_polysynapse",
        "has_electrical",
        "n_ex_test_spikes",
        "n_in_test_spikes",
        "distance",
        "lateral_distance",
        "vertical_distance",
    ], indent=4)

    if pair.has_synapse:
        print_heading(f"Synapse", level=2)
        syn = pair.synapse

        print_attrs(syn, [
            "synapse_type",
            "latency",
            "psc_amplitude",
            "psc_rise_time",
            "psc_decay_tau",
            "psp_amplitude",
            "psp_rise_time",
            "psp_decay_tau",
            "meta",
        ], indent=8)

        errs = [x for x in pipeline_jobs['synapse'].error.split('\n') if f'{pair_id[0]} {pair_id[1]}' in x]
        print("        pipeline errors:")
        for err in errs:
            print(" "*10 + err)





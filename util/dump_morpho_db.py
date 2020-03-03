import sys, json
from aisynphys.pipeline.morphology import import_morpho_db

target_file = sys.argv[1]

morpho_results = import_morpho_db()

json.dump(morpho_results, open(target_file, 'w'), indent=4)

import sys, json
from aisynphys import config

target_file = sys.argv[1]

cnxn_str = r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)}; DBQ=%s' % config.morpho_address
cnxn = pyodbc.connect(cnxn_str)
cursor = cnxn.cursor()
morpho_table = cursor.execute('select * from MPATCH_CellsofCluster')
results = morpho_table.fetchall()
morpho_results = {int(r.cell_specimen_id): r for r in results}

json.dump(morpho_results, open(target_file, 'w'))

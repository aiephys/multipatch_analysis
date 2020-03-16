import sys, json
import pyodbc
from aisynphys import config

target_file = sys.argv[1]

cnxn_str = r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)}; DBQ=%s' % config.morpho_address
cnxn = pyodbc.connect(cnxn_str)
cursor = cnxn.cursor()
morpho_table = cursor.execute('select * from MPATCH_CellsofCluster')
results = morpho_table.fetchall()
fields = [r[0] for r in results[0].cursor_description]
morpho_results = {r.cell_specimen_id: {k:getattr(r, k) for k in fields} for r in results}

json.dump(morpho_results, open(target_file, 'w'), indent=4)
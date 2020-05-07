# for making null lims_cells.json file in older experiments to submit to LIMS
import sys, os, json

directory = sys.argv[1]

json_file = os.path.join(directory, 'lims_cells.json')
if os.path.isfile(json_file):
    print('lims_cells.json file already exists in %s' % directory)
else:
    cells = []
    for i in range(1, 9):        
        cell = {
                'external_specimen_name': str(i),
                'patched_cell_container': None,
                'cell_reporter': None,
                'structure': None,
                }
        cells.append(cell)        

    json.dump(cells, open(json_file, 'w'))
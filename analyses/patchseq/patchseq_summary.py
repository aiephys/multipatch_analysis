from aisynphys.database import default_db as db
from aisynphys.connectivity import pair_was_probed
from aisynphys.ui.graphics import MatrixItem
import pyqtgraph as pg
import numpy as np
import time, random, sys, argparse

pg.dbg()

parser = argparse.ArgumentParser()
parser.add_argument('--timeline', action='store_true')
parser.add_argument('--hist', action='store_true')

args = parser.parse_args(sys.argv[1:])

species = {'mouse': ['o', (0, 0, 255, 150)],  'human': ['s', (0, 255, 0, 150)]}
metrics = ['area_400_10000bp', 'picogreen_yield','norm_marker_sum', 'genes_detected', 'tree_call']

q = db.query(db.PatchSeq, db.Experiment.acq_timestamp)
q = q.join(db.Cell, db.PatchSeq.cell_id==db.Cell.id)
q = q.join(db.Experiment, db.Cell.experiment_id==db.Experiment.id)
q = q.join(db.Slice, db.Experiment.slice_id==db.Slice.id)

for metric in metrics:
    ticks = None
    if args.timeline:
        time_plt = pg.plot()
        time_plt.addLegend()
    if args.hist:
        hist_plt = pg.plot()
        hist_plt.addLegend()
    for organism, symbols in species.items(): 
        patchseq_cells = q.filter(db.Slice.species==organism)
        symbol = symbols[0]
        color = symbols[1]

        values = [getattr(c[0], metric) for c in patchseq_cells if getattr(c[0], metric) is not None]
        if metric == 'tree_call':
            tree_map = {'Core': 4, 'I1': 3, 'I2': 2, 'I3': 1, 'PoorQ': 0}
            values = [tree_map[v] + random.uniform(-0.3, 0.3) for v in values]
            ticks = [[(tick_val, tick_name) for tick_name, tick_val in tree_map.items()]]
        if args.timeline and len(values) > 0:
            ts = [c[1] for c in patchseq_cells if getattr(c[0], metric) is not None]
            time_plt.plot(ts, values, symbol=symbol, symbolSize=10, pen=None, symbolBrush=color, name=organism)
            time_plt.setLabel('bottom', text='Experiment Date')
            time_plt.setLabel('left', metric)
            if ticks is not None:
                time_plt.getAxis('left').setTicks(ticks)
        if args.hist and len(values) > 0:
            if ticks is None:
                y, x = np.histogram(values, bins=np.linspace(0, max(values), 20))
            else:
                y, x = np.histogram(values, bins=len(list(tree_map.values())))
                hist_plt.getAxis('bottom').setTicks(ticks)
            hist_plt.plot(x, y, stepMode=True, fillLevel=0, brush=color, name=organism)
            hist_plt.setLabel('left', text='Number of Cells')
            hist_plt.setLabel('bottom', text=metric)
            mean_line = pg.InfiniteLine(pos=np.mean(values), angle=90, pen={'color': color[:3], 'width': 2})
            hist_plt.addItem(mean_line)

pre_q = db.query(db.Pair).join(db.PatchSeq, db.Pair.pre_cell_id==db.PatchSeq.cell_id).filter(db.PatchSeq.tree_call.in_(['Core', 'I1']))
pre_synapses = pre_q.filter(db.Pair.has_synapse==True).all()
print('%d synapses with well mapped presynaptic cells' % len(pre_synapses))

post_q = db.query(db.Pair).join(db.PatchSeq, db.Pair.post_cell_id==db.PatchSeq.cell_id).filter(db.PatchSeq.tree_call.in_(['Core', 'I1']))
post_synapses = post_q.filter(db.Pair.has_synapse==True).all()
print('%d synapses with well mapped postsynaptic cells' % len(post_synapses))

both_synapses = [pair for pair in pre_synapses if pair in post_synapses]
print('%d synapses with well mapped pre and postsynaptic cells' % len(both_synapses))

pre_pairs = pre_q.all()
post_pairs = post_q.all()
both_pairs = [pair for pair in pre_pairs if pair in post_pairs]

for organism in species.keys():
    print('%s Connectivity' % organism)
    cell_types = db.query(db.PatchSeq.tree_first_cluster).join(db.Cell).join(db.Experiment).join(db.Slice).filter(db.Slice.species==organism).filter(db.PatchSeq.tree_call.in_(['Core', 'I1'])).all()
    cell_types = set(cell_types)
    cell_types = sorted([ct[0] for ct in cell_types])
    pair_groups = {}
    for pre in cell_types:
        for post in cell_types:
            pair_groups[(pre, post)] = {'connected': None, 'probed': None}

    for pair in both_pairs:
        if pair.experiment.slice.species != organism:
            continue
        pre_type = pair.pre_cell.patch_seq.tree_first_cluster
        post_type = pair.post_cell.patch_seq.tree_first_cluster
        if organism == 'mouse':
            synapse_type = 'ex' if pre_type.startswith('L2/3') else 'in'
        if organism == 'human':
            synapse_type = 'ex' if pre_type.startswith('EXC') else 'in'
        pair_probed = pair_was_probed(pair, synapse_type)
        if pair_probed:
            pair_groups[(pre_type, post_type)]['probed'] = 1 if pair_groups[(pre_type, post_type)]['probed'] is None else pair_groups[(pre_type, post_type)]['probed'] + 1
            if pair_groups[(pre_type, post_type)]['connected'] is None:
                pair_groups[(pre_type, post_type)]['connected'] = 0
            if pair.has_synapse:
                pair_groups[(pre_type, post_type)]['connected'] += 1

    if organism ==  'mouse':
        shape = (len(cell_types), len(cell_types))
        text = np.empty(shape, dtype=object)
        text.fill('')
        fgcolor = np.empty(shape, dtype=object)
        fgcolor.fill(0.6)
        bgcolor = np.empty(shape, dtype=object)
        bgcolor.fill(tuple(np.array([128,128,128,255])))
        bordercolor = np.empty(shape, dtype=object)
        bordercolor.fill(0.8)
        cmap = pg.ColorMap(
                [0, 0.01, 0.03, 0.1, 0.3, 1.0],
                [(0,0,100, 255), (80,0,80, 255), (140,0,0, 255), (255,100,0, 255), (255,255,100, 255), (255,255,255, 255)]
                )
        
        for i, pre in enumerate(cell_types):
            for j, post in enumerate(cell_types):
                data = pair_groups[(pre, post)]
                probed = data['probed']
                if probed is None:
                    continue
                conn = data['connected']
                cprob = conn/probed
                text[i, j] = "%d/%d" % (conn, probed)
                color = cmap.map(cprob)
                bgcolor[i, j] = color
                fgcolor[i, j] = 'w' if sum(color) < 384 else 'k'

        matrix = MatrixItem(text=text, fgcolor=fgcolor, bgcolor=bgcolor, rows=cell_types, cols=cell_types, border_color=bordercolor, header_color='k')
        window = pg.GraphicsLayoutWidget()
        viewbox = window.addViewBox()
        viewbox.addItem(matrix)
        viewbox.invertY()
        viewbox.setAspectLocked()
        viewbox.setBackgroundColor('w')
        window.show()

    total_conn = sum([group['connected'] for group in pair_groups.values() if group['probed'] is not None])
    total_probed = sum([group['probed'] for group in pair_groups.values() if group['probed'] is not None])
    print('Total Connected / Total Probed:\t%d/%d' % (total_conn, total_probed))


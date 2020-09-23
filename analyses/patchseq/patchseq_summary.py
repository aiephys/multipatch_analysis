from aisynphys.database import default_db as db
from aisynphys.connectivity import pair_was_probed
from aisynphys.ui.graphics import MatrixItem
from datetime import datetime
import pyqtgraph as pg
import numpy as np
import time, random, sys, argparse

pg.dbg()

parser = argparse.ArgumentParser()
parser.add_argument('--timeline', action='store_true')
parser.add_argument('--hist', action='store_true')
parser.add_argument('--date-range', nargs=2)
parser.add_argument('--species')

args = parser.parse_args(sys.argv[1:])

species = {'mouse': ['o', (0, 0, 255, 150)],  'human': ['s', (0, 255, 0, 150)]}
if args.species is not None:
    species = {args.species: species[args.species]} 
metrics = ['area_400_10000bp', 'picogreen_yield','norm_marker_sum', 'genes_detected', 'tree_call']
subclass_colors = {'pvalb': (0, 0, 255, 150), 'sst': (255, 0, 0, 150), 'vip': (0, 255, 0, 150), 'unknown': (150, 150, 150, 150), 'other': (150, 0, 150, 150)}

q = db.query(db.PatchSeq, db.Experiment.acq_timestamp)
q = q.join(db.Cell, db.PatchSeq.cell_id==db.Cell.id)
q = q.join(db.Experiment, db.Cell.experiment_id==db.Experiment.id)
q = q.join(db.Slice, db.Experiment.slice_id==db.Slice.id)

if args.date_range is not None:
    dates = args.date_range
    start = datetime.strptime(dates[0], "%Y/%m/%d")
    end = datetime.strptime(dates[1], "%Y/%m/%d")
    q = q.filter(db.Experiment.date >= start)
    q = q.filter(db.Experiment.date <= end)

for metric in metrics:
    ticks = None
    no_data = False
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
        labels = None

        values = [getattr(c[0], metric) for c in patchseq_cells if getattr(c[0], metric) is not None]
        
        if len(values)==0:
            no_data = True
            continue
        if len(species)==1 and organism=='mouse':
            color = [pg.mkColor(subclass_colors.get(c[0].cell.cre_type, subclass_colors['other'])) for c in patchseq_cells if getattr(c[0], metric) is not None]
            labels = [pg.TextItem(text=key, color=value[:3]) for key, value in subclass_colors.items()]
        no_data = False
        if metric == 'tree_call':
            tree_map = {'Core': 4, 'I1': 3, 'I2': 2, 'I3': 1, 'PoorQ': 0}
            values = [tree_map[v] + random.uniform(-0.3, 0.3) for v in values]
            ticks = [[(tick_val, tick_name) for tick_name, tick_val in tree_map.items()]]
        if args.timeline and len(values) > 0:
            ts = [c[1] for c in patchseq_cells if getattr(c[0], metric) is not None]
            days = [datetime.fromtimestamp(ts_val).date() for ts_val in ts]
            unique_days = [ts[days.index(x)] for x in set(days)]
            date_ticks = [[(ts_val, datetime.fromtimestamp(ts_val).strftime("%m\n%d\n%y")) for ts_val in unique_days]]
            time_plt.plot(ts, values, symbol=symbol, symbolSize=10, pen=None, symbolBrush=color, name=organism)
            time_plt.setLabel('bottom', text='Experiment Date')
            time_plt.setLabel('left', metric)
            time_plt.getAxis('bottom').setTicks(date_ticks)
            if labels is not None:
                label_items = [time_plt.addItem(label) for label in labels]
                spacing = np.mean(np.diff(values)) *50 
                [label.setPos(max(ts), max(values)-i*spacing) for i,  label in enumerate(labels)] 
            if ticks is not None:
                time_plt.getAxis('left').setTicks(ticks)
        if args.hist and len(values) > 0:
            if ticks is None:
                y, x = np.histogram(values, bins=np.linspace(0, max(values), 20))
            else:
                y, x = np.histogram(values, bins=len(list(tree_map.values())))
                hist_plt.getAxis('bottom').setTicks(ticks)
            hist_plt.plot(x, y, stepMode=True, fillLevel=0, brush=symbols[1], name=organism)
            hist_plt.setLabel('left', text='Number of Cells')
            hist_plt.setLabel('bottom', text=metric)
            mean_line = pg.InfiniteLine(pos=np.nanmean(values), angle=90, pen={'color': symbols[1], 'width': 2})
            hist_plt.addItem(mean_line)
    if no_data == True:
        try:
            time_plt.win.close()
        except NameError:
            continue

        try:
            hist_plt.win.close()
        except NameError:
            continue

pair_q = db.query(db.Pair)
if args.date_range is not None:
    dates = args.date_range
    start = datetime.strptime(dates[0], "%Y/%m/%d")
    end = datetime.strptime(dates[1], "%Y/%m/%d")
    pair_q = pair_q.join(db.Experiment, db.Pair.experiment_id==db.Experiment.id).filter(db.Experiment.date >= start)
    pair_q = pair_q.filter(db.Experiment.date <= end)

pre_q = pair_q.join(db.PatchSeq, db.Pair.pre_cell_id==db.PatchSeq.cell_id).filter(db.PatchSeq.tree_call.in_(['Core', 'I1']))
pre_synapses = pre_q.filter(db.Pair.has_synapse==True).all()
print('%d synapses with well mapped presynaptic cells' % len(pre_synapses))

post_q = pair_q.join(db.PatchSeq, db.Pair.post_cell_id==db.PatchSeq.cell_id).filter(db.PatchSeq.tree_call.in_(['Core', 'I1']))
post_synapses = post_q.filter(db.Pair.has_synapse==True).all()
print('%d synapses with well mapped postsynaptic cells' % len(post_synapses))

both_synapses = [pair for pair in pre_synapses if pair in post_synapses]
print('%d synapses with well mapped pre and postsynaptic cells' % len(both_synapses))

cmap = pg.ColorMap(
        [0, 0.01, 0.03, 0.1, 0.3, 1.0],
        [(0,0,100, 255), (80,0,80, 255), (140,0,0, 255), (255,100,0, 255), (255,255,100, 255), (255,255,255, 255)]
        )

pre_pairs = pre_q.all()
post_pairs = post_q.all()
both_pairs = [pair for pair in pre_pairs if pair in post_pairs]


for organism in species.keys():
    print('%s Connectivity' % organism)
    cells = q.filter(db.Slice.species==organism).filter(db.PatchSeq.tree_call.in_(['Core', 'I1'])).all()
    cell_types = [c[0].tree_first_cluster for c in cells]
    cell_types = set(cell_types)
    # cell_types = sorted([ct[0] for ct in cell_types])
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

    total_conn = sum([group['connected'] for group in pair_groups.values() if group['probed'] is not None])
    total_probed = sum([group['probed'] for group in pair_groups.values() if group['probed'] is not None])
    print('Total Connected / Total Probed:\t%d/%d' % (total_conn, total_probed))

cell_types = db.query(db.PatchSeq.tree_first_cluster).join(db.Cell).join(db.Experiment).join(db.Slice).filter(db.Slice.species=='mouse').filter(db.PatchSeq.tree_call.in_(['Core', 'I1'])).all()
cell_types = set(cell_types)
cell_types = sorted([ct[0] for ct in cell_types])

## mouse t-types matrix

pair_groups = {}
for pre in cell_types:
    for post in cell_types:
        pair_groups[(pre, post)] = {'connected': None, 'probed': None}

for pair in both_pairs:
    if pair.experiment.slice.species != 'mouse':
        continue
    synapse_type = 'ex' if pre_type.startswith('L2/3') else 'in'
    pair_probed = pair_was_probed(pair, synapse_type)
    pre_type = pair.pre_cell.patch_seq.tree_first_cluster
    post_type = pair.post_cell.patch_seq.tree_first_cluster
    if pair_probed:
        pair_groups[(pre_type, post_type)]['probed'] = 1 if pair_groups[(pre_type, post_type)]['probed'] is None else pair_groups[(pre_type, post_type)]['probed'] + 1
        if pair_groups[(pre_type, post_type)]['connected'] is None:
            pair_groups[(pre_type, post_type)]['connected'] = 0
        if pair.has_synapse:
            pair_groups[(pre_type, post_type)]['connected'] += 1

shape = (len(cell_types), len(cell_types))
text = np.empty(shape, dtype=object)
text.fill('')
fgcolor = np.empty(shape, dtype=object)
fgcolor.fill(0.6)
bgcolor = np.empty(shape, dtype=object)
bgcolor.fill(tuple(np.array([128,128,128,255])))
bordercolor = np.empty(shape, dtype=object)
bordercolor.fill(0.8)


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

if np.all(fgcolor==0.6) == False:
    matrix = MatrixItem(text=text, fgcolor=fgcolor, bgcolor=bgcolor, rows=cell_types, cols=cell_types, border_color=bordercolor, header_color='k')
    window = pg.GraphicsLayoutWidget()
    window.show()
    viewbox = window.addViewBox()
    viewbox.addItem(matrix)
    viewbox.invertY()
    viewbox.setAspectLocked()
    viewbox.setBackgroundColor('w')


## score summary
# tree_scores = ['tree_first_score', 'tree_second_score', 'tree_first_kl', 'tree_second_kl', 'tree_first_cor', 'tree_second_cor']

# cells = q.filter(db.PatchSeq.tree_call != None).all()
# tree_map = {'Core': 4, 'I1': 3, 'I2': 2, 'I3': 1, 'PoorQ': 0}
# ticks = [[(tick_val, tick_name) for tick_name, tick_val in tree_map.items()]]

# for name in tree_scores:
#     scores = [getattr(cell[0], name) for cell in cells if hasattr(cell[0], name)]
#     tree_call = [cell[0].tree_call for cell in cells if hasattr(cell[0], name)]
#     tree_num = [tree_map[v]+random.uniform(-0.3, 0.3) for v in tree_call]
#     plt = pg.plot(tree_num, scores, symbol='o', pen=None, symbolBrush=(0, 0, 255, 150))
#     plt.getAxis('bottom').setTicks(ticks)
#     if name.endswith('score'):
#         name = '_'.join(name.split('_')[0:2]+['bt'])
#     plt.setLabel('left', text=name)

plt = pg.plot()
no_data = False
for organism, symbols in species.items():
    symbol = symbols[0]
    color = symbols[1]
    name = organism
    cells = q.filter(db.PatchSeq.tree_call != None).filter(db.Slice.species==organism).all()
    if len(cells) == 0:
        no_data = True
        continue
    tree_call = [cell[0].tree_call for cell in cells]
    if len(species)==1 and organism=='mouse':
        color = [pg.mkColor(subclass_colors.get(c[0].cell.cre_type, subclass_colors['other'])) for c in cells]
        name = subclass_colors.keys()
    tree_map = {'Core': 4, 'I1': 3, 'I2': 2, 'I3': 1, 'PoorQ': 0}
    nucleus = [cell[0].nucleus for cell in cells]
    nuc_map = {True: 1, False: 0, None: -1}
    tree_num = [tree_map[v]+random.uniform(-0.3, 0.3) for v in tree_call]
    if organism == 'mouse':     
        nuc_num = [nuc_map[n] + random.uniform(0, 0.3) for n in nucleus]
    else:
        nuc_num = [nuc_map[n] + random.uniform(-0.3, 0) for n in nucleus]
    no_data = False
    plt.plot(tree_num, nuc_num, pen=None, symbol=symbol, symbolBrush=color, name=name)
plt.getAxis('bottom').setTicks(ticks)
plt.getAxis('left').setTicks([[(-1, 'No Data'), (0, 'Nucleus Absent'), (1, 'Nucleus Present')]])
if no_data == True:
    plt.win.close()

cells = q.all()
expts = [c[1] for c in cells]
expts = db.query(db.Experiment).filter(db.Experiment.acq_timestamp.in_(expts)).all()
patched=[]
nucleus=[]
for e in expts:
    patched.append(len(e.cells.keys()))
    nuc =0
    for cell in e.cells.values():
        if cell.patch_seq is not None:
            if cell.patch_seq.nucleus is True:
                nuc+=1
    nucleus.append(nuc)

p = np.array(patched)
n = np.array(nucleus)
n1 = n/p

p_scatter = [p1+random.uniform(-0.2, 0.2) for p1 in p]
# n1_scatter = [n2+random.uniform(-0.2, 0.2) for n2 in n1]

plt = pg.plot(p_scatter, n1, pen=None, symbol='o')
plt.setLabel('bottom', text='Number of cells patched')
plt.setLabel('left', text='Fraction of nuclei recovered')

## human subclasses

if 'human' in species.keys():
    subclass = db.query(db.PatchSeq.subclass_label).join(db.Cell).join(db.Experiment).join(db.Slice).filter(db.Slice.species=='human').filter(db.PatchSeq.norm_marker_sum>=0.4).all()
    subclass = set(subclass)
    subclass = sorted([ct[0] for ct in subclass])

    pre_q = pair_q.join(db.PatchSeq, db.Pair.pre_cell_id==db.PatchSeq.cell_id).filter(db.PatchSeq.norm_marker_sum>=0.4)
    post_q = pair_q.join(db.PatchSeq, db.Pair.post_cell_id==db.PatchSeq.cell_id).filter(db.PatchSeq.norm_marker_sum>=0.4)

    pre_pairs = pre_q.all()
    post_pairs = post_q.all()
    both_pairs = [pair for pair in pre_pairs if pair in post_pairs]

    window2 = pg.GraphicsLayoutWidget()
    window2.show()

    pair_groups = {}
    for pre in subclass:
        for post in subclass:
            pair_groups[(pre, post)] = {'connected': None, 'probed': None}

    for pair in both_pairs:
        if pair.experiment.slice.species != 'human':
            continue
        synapse_type = 'in' if pre_type in ['SST', 'PVALB', 'VIP', 'LAMP5/PAX6/Other'] else 'ex'
        pair_probed = pair_was_probed(pair, synapse_type)
        pre_type = pair.pre_cell.patch_seq.subclass_label
        post_type = pair.post_cell.patch_seq.subclass_label
        if pair_probed:
            pair_groups[(pre_type, post_type)]['probed'] = 1 if pair_groups[(pre_type, post_type)]['probed'] is None else pair_groups[(pre_type, post_type)]['probed'] + 1
            if pair_groups[(pre_type, post_type)]['connected'] is None:
                pair_groups[(pre_type, post_type)]['connected'] = 0
            if pair.has_synapse:
                pair_groups[(pre_type, post_type)]['connected'] += 1

    shape = (len(subclass), len(subclass))
    text = np.empty(shape, dtype=object)
    text.fill('')
    fgcolor = np.empty(shape, dtype=object)
    fgcolor.fill(0.6)
    bgcolor = np.empty(shape, dtype=object)
    bgcolor.fill(tuple(np.array([128,128,128,255])))
    bordercolor = np.empty(shape, dtype=object)
    bordercolor.fill(0.8)

            
    for i, pre in enumerate(subclass):
        for j, post in enumerate(subclass):
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

    matrix = MatrixItem(text=text, fgcolor=fgcolor, bgcolor=bgcolor, rows=subclass, cols=subclass, border_color=bordercolor, header_color='k')

    viewbox = window2.addViewBox()
    viewbox.addItem(matrix)
    viewbox.invertY()
    viewbox.setAspectLocked()
    viewbox.setBackgroundColor('w')
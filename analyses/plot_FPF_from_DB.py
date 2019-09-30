import argparse, first_pulse_from_DB
import pyqtgraph as pg
from sqlalchemy.orm import aliased
from neuroanalysis.ui.plot_grid import PlotGrid
from aisynphys.database import database as db
from neuroanalysis.data import TSeries, TSeriesList

def plot_features(organism=None, conn_type=None, calcium=None, age=None, sweep_thresh=None, fit_thresh=None):
    s = db.session()

    filters = {
        'organism': organism,
        'conn_type': conn_type,
        'calcium': calcium,
        'age': age
    }

    selection = [{}]
    for key, value in filters.iteritems():
        if value is not None:
            temp_list = []
            value_list = value.split(',')
            for v in value_list:
                temp = [s1.copy() for s1 in selection]
                for t in temp:
                    t[key] = v
                temp_list = temp_list + temp
            selection = list(temp_list)

    if len(selection) > 0:

        response_grid = PlotGrid()
        response_grid.set_shape(len(selection), 1)
        response_grid.show()
        feature_grid = PlotGrid()
        feature_grid.set_shape(6, 1)
        feature_grid.show()

        for i, select in enumerate(selection):
            pre_cell = aliased(db.Cell)
            post_cell = aliased(db.Cell)
            q_filter = []
            if sweep_thresh is not None:
                q_filter.append(FirstPulseFeatures.n_sweeps>=sweep_thresh)
            species = select.get('organism')
            if species is not None:
                q_filter.append(db.Slice.species==species)
            c_type = select.get('conn_type')
            if c_type is not None:
                pre_type, post_type = c_type.split('-')
                pre_layer, pre_cre = pre_type.split(';')
                if pre_layer == 'None':
                    pre_layer = None
                post_layer, post_cre = post_type.split(';')
                if post_layer == 'None':
                    post_layer = None
                q_filter.extend([pre_cell.cre_type==pre_cre, pre_cell.target_layer==pre_layer,
                                post_cell.cre_type==post_cre, post_cell.target_layer==post_layer])
            calc_conc = select.get('calcium')
            if calc_conc is not None:
                q_filter.append(db.Experiment.acsf.like(calc_conc + '%'))
            age_range = select.get('age')
            if age_range is not None:
                age_lower, age_upper = age_range.split('-')
                q_filter.append(db.Slice.age.between(int(age_lower), int(age_upper)))

            q = s.query(FirstPulseFeatures).join(db.Pair, FirstPulseFeatures.pair_id==db.Pair.id)\
                .join(pre_cell, db.Pair.pre_cell_id==pre_cell.id)\
                .join(post_cell, db.Pair.post_cell_id==post_cell.id)\
                .join(db.Experiment, db.Experiment.id==db.Pair.expt_id)\
                .join(db.Slice, db.Slice.id==db.Experiment.slice_id)

            for filter_arg in q_filter:
                q = q.filter(filter_arg)

            results = q.all()

            trace_list = []
            for pair in results:
                #TODO set t0 to latency to align to foot of PSP
                trace = TSeries(data=pair.avg_psp, sample_rate=db.default_sample_rate)
                trace_list.append(trace)
                response_grid[i, 0].plot(trace.time_values, trace.data)
            if len(trace_list) > 0:
                grand_trace = TSeriesList(trace_list).mean()
                response_grid[i, 0].plot(grand_trace.time_values, grand_trace.data, pen='b')
                response_grid[i, 0].setTitle('layer %s, %s-> layer %s, %s; n_synapses = %d' %
                                          (pre_layer, pre_cre, post_layer, post_cre, len(trace_list)))
            else:
                print('No synapses for layer %s, %s-> layer %s, %s' % (pre_layer, pre_cre, post_layer, post_cre))

    return response_grid, feature_grid

if __name__ == '__main__':
    selection = parser.add_argument_group('Selection criteria', 'Filters for experiment analysis')
    selection.add_argument('--organism', type=str, help='mouse, human, or mouse,human will show both')
    selection.add_argument('--conn-type', type=str, help="""Enter as layer;cre-pre-layer;cre-post.
                                                      If using multiple connection types separate with ",".
                                                      Ex 5;pvalb-5;pavalb,5;pvalb-5;sst""")
    selection.add_argument('--calcium', type=str , help='Enter calcium concentration, may enter multiple separated by ","')
    selection.add_argument('--age', type=str, help="""Enter the age range separate by "-". May enter multiple ranges
                                                separated by ",". Will only apply to mouse """)
    selection.add_argument('--sweep-thresh', type=int, help='Minimum number of sweeps required for analysis')
    selection.add_argument('--fit-thresh', type=float, help="Discard results who's NRMSE from the fit exeed this value")

    response_plot, feature_plot = plot_features(organism=args.organism, conn_type=args.conn_type, calcium=args.calcium,
                                                age=args.age,
                                                sweep_thresh=args.sweep_thresh, fit_thresh=args.fit_thresh)



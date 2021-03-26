"""
Analysis code used to load tables of synapse features and perform dimensionality reduction on them,
with the intent to better understand what the major classes of synapses are, and how they relate
to cell and synapse properties.
"""
import os, pickle
import numpy as np
import pandas
import sklearn.preprocessing, sklearn.pipeline
import umap
from .database import default_db
from neuroanalysis.util.optional_import import optional_import
plt = optional_import('matplotlib.pyplot')
sns = optional_import('seaborn')


# Complete list of fields to load along with synapses by default
default_data_fields = [
    'experiment_ext_id',
    'pre_ext_id',
    'post_ext_id',
    'synapse_type',
    'pre_cre_type', 
    'post_cre_type', 
    'pre_cell_class', 
    'post_cell_class', 
    'species',
    'pre_layer', 
    'post_layer',
    'distance',
    'n_ex_test_spikes',
    'n_in_test_spikes',
    'psp_amplitude', 
    'psc_amplitude', 
    'psp_rise_time', 
    'psp_decay_tau', 
    'psc_rise_time', 
    'psc_decay_tau',
    'latency',
    'paired_pulse_ratio_50hz',
    'stp_initial_50hz',
    'stp_initial_50hz_n',
    'stp_initial_50hz_std',
    'stp_induction_50hz',
    'stp_induction_50hz_n',
    'stp_induction_50hz_std',
    'stp_recovery_250ms',
    'stp_recovery_250ms_n',
    'stp_recovery_250ms_std',
    'stp_recovery_single_250ms',
    'stp_recovery_single_250ms_n',
    'pulse_amp_90th_percentile',
    'noise_amp_90th_percentile',
    'noise_std',
    'variability_resting_state',
    'variability_second_pulse_50hz',
    'variability_stp_induced_state_50hz',
    'variability_change_initial_50hz',
    'variability_change_induction_50hz',
    'paired_event_correlation_1_2_r',
    'paired_event_correlation_1_2_p',
    'paired_event_correlation_2_4_r',
    'paired_event_correlation_2_4_p',
    'paired_event_correlation_4_8_r',
    'paired_event_correlation_4_8_p',
    'n_model_source_events',
]



def synapse_query(db=default_db):
    """Return a query that selects synapses joined to various cell and synapse properties.
    """
    q = db.pair_query(
        synapse=True, 
        project_name=['mouse V1 coarse matrix', 'mouse V1 pre-production', 'human coarse matrix'], 
    )
    pre_cell = q.pre_cell
    post_cell = q.post_cell
    q = (q
        .join(db.SynapseModel)
        .add_entity(db.Synapse)
        .add_entity(db.Dynamics)
        .add_column(pre_cell.cre_type.label('pre_cre_type'))
        .add_column(post_cell.cre_type.label('post_cre_type'))
        .add_column(pre_cell.cell_class.label('pre_cell_class'))
        .add_column(post_cell.cell_class.label('post_cell_class'))
        .add_column(db.Slice.species)
        .add_column(db.Experiment.ext_id.label('experiment_ext_id'))
        .add_column(pre_cell.ext_id.label('pre_ext_id'))
        .add_column(post_cell.ext_id.label('post_ext_id'))
        .add_column(q.pre_location.cortical_layer.label('pre_layer'))
        .add_column(q.post_location.cortical_layer.label('post_layer'))
        .add_column(db.SynapseModel.n_source_events.label('n_model_source_events'))
    )
    return q


def load_data(query=None, data_fields=default_data_fields):
    if query is None:
        query = synapse_query()

    synapses = query.dataframe()
    syn_data = synapses[data_fields].copy()

    return syn_data.set_index(['experiment_ext_id', 'pre_ext_id', 'post_ext_id'])


def label_subclasses(data, subclasses):
    """Return a copy of *data* with pre_ and post_subclass columns added based on the definitions provided
    in *subclasses*.

    Example::

        subclasses = {
            'Pvalb': {'cre_type': 'pvalb'},
            'Sst': {'cre_type': 'sst'},
            'Vip': {'cre_type': 'vip'},
            'L2/3E': {'species': 'mouse', 'layer': '2/3', 'cell_class': 'ex', 'cre_type': 'unknown'},
        }            
    """
    data = data.copy()

    # label rows with subclass identifiers
    data['pre_subclass'] = ""
    data['post_subclass'] = ""

    # keys that don't have pre_ or post_ prepended
    pass_keys = ['species']

    for side in 'pre_', 'post_':
        sided_classes = {}
        # prepend some keys with pre_ or post_
        for name,crit in subclasses.items():
            sided_classes[name] = {(k if k in pass_keys else side+k):v  for k, v in crit.items()}
        # label rows based on pre or post criteria
        data = label_rows(data, sided_classes, side+'subclass')

    return data


def label_rows(data, subclasses, new_column):
    """Return a copy of *data* with a *new_column* added that labels
    each row based on the criteria listed in *subclasses*.
    """
    data2 = data.copy()
    for cls_name, cls_crit in subclasses.items():
        mask = np.ones(len(data), dtype=bool)
        for k,vals in cls_crit.items():
            if not isinstance(vals, list):
                vals = [vals]
            mask2 = np.zeros(len(data), dtype=bool)
            for v in vals:
                mask2 |= data[k] == v        
            mask &= mask2
        data2.loc[mask, new_column] = cls_name
    return data2


def load_model_vectors(model_file, max_vector_size=np.inf):
    """Return a pandas dataframe containing sPCA vectors describing the posterior distribution of likelihood
    values across model parameters.

    The *model_file* is generated by aisynphys/analyses/stochastic_model_reduction.py
    """
    sm_results = pickle.load(open(model_file, 'rb'))
    n_cols = min(max_vector_size, sm_results['result'].shape[1])

    results = []
    for i,cache_file in enumerate(sm_results['cache_files']):
        expt_id, pre_cell_id, post_cell_id = os.path.split(os.path.splitext(cache_file)[0])[1].split('_')
        row = {
            'experiment_ext_id': expt_id,
            'pre_ext_id': pre_cell_id,
            'post_ext_id': post_cell_id,
        }
        for j in range(n_cols):
            row['PC_%d'%j] = sm_results['result'][i, j]
        results.append(row)

    return pandas.DataFrame(results).set_index(['experiment_ext_id', 'pre_ext_id', 'post_ext_id'])


def umap_pipeline(n_components=2, n_neighbors=15, min_dist=0.8, random_state=0, **kwds):
    """Create an sklearn Pipeline containing a power transfomer (for normalization) and UMAP.
    All arguments are used to construct the UMAP instance.
    """
    normalizer = sklearn.preprocessing.PowerTransformer(method='yeo-johnson', standardize=True)
    mapper = sklearn.pipeline.Pipeline([
        ('normalize', normalizer),
        ('reduce', umap.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            random_state=random_state,
            **kwds
        )),
    ])
    return mapper


def run_umap_pipeline(data, mapper, features):
    """Filter *data* to include only rows that contain values for all *features*, 
    then train *mapping* on this filtered data and return a a dataframe containing
    the filtered data plus umap embedding columns.
    """

    # only keep synapses with all feature fields
    mask = data[features].isna().any(axis=1)

    clean_syn_data = data.loc[~mask].copy()

    print("Dropped to %d synapses" % len(clean_syn_data))

    feature_data = clean_syn_data[features]
    mapper.fit(feature_data)

    embedding = pandas.DataFrame(
        data=mapper.transform(feature_data),
        index=clean_syn_data.index,
        columns=['umap-0', 'umap-1'],
    )

    return clean_syn_data.join(embedding)


def show_umap(data, ax, title=None, legend_title=None, picking=True, **kwds):
    """Default styling and convenience features for generating umap scatter plots.
    """
    x, y = 'umap-0', 'umap-1'

    opts = dict(
        x=x, y=y,
        linewidth=0,
        s=70,
        ax=ax,
        legend=False,
        picker=True,
    )
    opts.update(kwds)
    
    # first show all points in grey
    sns.scatterplot(data=data, x=x, y=y, color=(0, 0, 0, 0.1), linewidth=0, s=opts['s'] / 3, ax=ax)

    if 'palette' in opts and isinstance(opts['palette'], dict):
        # automatically order legend when we provide a dictionary palette
        if 'hue_order' not in opts:
            opts['hue_order'] = opts['palette'].keys()
        # and mask out data that is missing an entry in palette
        mask = np.array([x in opts['palette'] for x in data[opts['hue']]])        
        data = data[mask]

    sp = sns.scatterplot(data=data, **opts)
    sns.despine(ax=ax, top=True, right=True, left=True, bottom=True, offset=None, trim=False)
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)

    if title is not None:
        ax.set_title(title)
    if legend_title is not None:
        sp.legend().set_title(legend_title)

    if picking:
        text = ax.text(0, 1, "", va="top", ha="left", transform=ax.transAxes)

        def onpick(event):
            t = '\n'.join([" ".join(data.index[i]) for i in event.ind])
            text.set_text(t)
            print(t)

        cid = ax.figure.canvas.mpl_connect('pick_event', onpick)

    return sp


def show_subclass_umap(data, subclasses, palette, ax, **kwds):
    data2 = label_rows(data, subclasses, new_column='subclass_column')
    show_umap(data2, hue="subclass_column", palette=palette, ax=ax, **kwds)



from multipatch_analysis.database import database as db
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.fitting import fit_psp
from neuroanalysis.baseline import float_mode

import argparse, sys, strength_analysis
import pyqtgraph as pg
import numpy as np


class TableGroup(object):
    def __init__(self):
        self.mappings = {}
        self.create_mappings()

    def __getitem__(self, item):
        return self.mappings[item]

    def create_mappings(self):
        for k,schema in self.schemas.items():
            self.mappings[k] = db.generate_mapping(k, schema)

    def drop_tables(self):
        for k in self.schemas:
            if k in db.engine.table_names():
                self[k].__table__.drop(bind=db.engine)

    def create_tables(self):
        for k in self.schemas:
            if k not in db.engine.table_names():
                self[k].__table__.create(bind=db.engine)


class FirstPulseFeaturesTableGroup(TableGroup):
    schemas = {
        'first_pulse_features': [
            """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed qc in current clamp""",
            ('pair_id', 'pair.id', '', {'index': True}),
            ('ic_fit_amp', 'float', 'amplitude from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_latency', 'float', 'latency from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_rise_time', 'float', 'rise time from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_decay_tau', 'float', 'decay tau from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_amp_cv', 'float', 'coefficient of variation for first pulse amplitude in 10, 20, 50 Hz stimuli'),
            ('avg_psp', 'array', 'average psp time series, spike aligned, baseline subtracted'),
            ('n_sweeps', 'int', 'number of sweeps in avg_psp'),
            ('pulse_ids', 'object', 'list of first pulse ids in avg_psp, len(pulse_ids) == n_sweeps')
            #('ic_fit_NRMSE', 'float', 'NRMSE returned from psp_fit')
            #TODO: consider removing 50Hz responses from decay calculation
        ]
    }

    def create_mappings(self):
        TableGroup.create_mappings(self)

        FirstPulseFeatures = self['first_pulse_features']

        db.Pair.first_pulse_features = db.relationship(FirstPulseFeatures, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        FirstPulseFeatures.pair = db.relationship(db.Pair, back_populates="first_pulse_features", single_parent=True)

first_pulse_features_tables = FirstPulseFeaturesTableGroup()


def init_tables():
    global FirstPulseFeatures
    first_pulse_features_tables.create_tables()
    FirstPulseFeatures = first_pulse_features_tables['first_pulse_features']

def update_analysis(limit=None):
    s = db.Session()
    q = s.query(db.Pair, FirstPulseFeatures).outerjoin(FirstPulseFeatures).filter(FirstPulseFeatures.pair_id == None)
    if limit is not None:
        q = q.limit(limit)
    print("Updating %d pairs.." % q.count())
    records = q.all()
    for i,record in enumerate(records):
        pair = record[0]
        pulse_responses, pulse_ids, pulse_response_amps = filter_pulse_responses(pair)
        if len(pulse_responses) > 0:
            results = first_pulse_features(pair, pulse_responses, pulse_response_amps)
            fpf = FirstPulseFeatures(pair=pair, n_sweeps=len(pulse_ids), pulse_ids=pulse_ids, **results)
            s.add(fpf)
            if i % 10 == 0:
                s.commit()
                print("%d pairs added to the DB of %d" %(i, len(records)))
    s.commit()
    s.close()


def first_pulse_features(pair, pulse_responses, pulse_response_amps):

    avg_psp = TraceList(pulse_responses).mean()
    dt = avg_psp.dt
    avg_psp_baseline = float_mode(avg_psp.data[:int(10e-3/dt)])
    avg_psp_bsub = avg_psp.copy(data=avg_psp.data - avg_psp_baseline)
    lower_bound = -float('inf')
    upper_bound = float('inf')
    xoffset = pair.connection_strength.ic_fit_xoffset
    if xoffset is None:
        xoffset = 14*10e-3
    synapse_type = pair.connection_strength.synapse_type
    if synapse_type == 'ex':
        amp_sign = '+'
    elif synapse_type == 'in':
        amp_sign = '-'
    else:
        raise Exception('Synapse type is not defined, reconsider fitting this pair %s %d->%d' %
                        (pair.expt_id, pair.pre_cell_id, pair.post_cell_id))

    weight = np.ones(len(avg_psp.data)) * 10.  # set everything to ten initially
    weight[int(10e-3 / dt):int(12e-3 / dt)] = 0.  # area around stim artifact
    weight[int(12e-3 / dt):int(19e-3 / dt)] = 30.  # area around steep PSP rise

    psp_fits = fit_psp(avg_psp,
                       xoffset=(xoffset, lower_bound, upper_bound),
                       yoffset=(avg_psp_baseline, lower_bound, upper_bound),
                       sign=amp_sign,
                       weight=weight)

    amp_cv = np.std(pulse_response_amps)/np.mean(pulse_response_amps)

    features = {'ic_fit_amp': psp_fits.best_values['amp'],
                'ic_fit_latency': psp_fits.best_values['xoffset'] - 10e-3,
                'ic_fit_rise_time': psp_fits.best_values['rise_time'],
                'ic_fit_decay_tau': psp_fits.best_values['decay_tau'],
                'ic_amp_cv': amp_cv,
                'avg_psp': avg_psp_bsub.data}
                #'ic_fit_NRMSE': psp_fits.nrmse()} TODO: nrmse not returned from psp_fits?

    return features

def filter_pulse_responses(pair):
    ### get first pulse response if it passes qc for excitatory or inhibitory analysis

    # TODO: learn how to do what's below in one query
    # s = db.Session()
    # q = s.query(db.PulseResponse.data, db.StimSpike, db.PatchClampRecording)
    # q = q.join(db.StimPulse).join(db.StimSpike).join(db.PatchClampRecording)
    # filters = [
    #     (db.Pair == pair)
    #     (db.StimPulse.pulse_number == 1),
    #     (db.StimPulse.n_spikes == 1),
    #     (db.StimSpike.max_dvdt_time != None),
    #     (db.PulseResponse.ex_qc_pass == True)
    #     (db.PatchClampRecording.clamp_mode == 'ic')
    # ]
    #
    # for filter_arg in filters:
    #     q = q.filter(*filter_arg)

    synapse_type = pair.connection_strength.synapse_type
    pulse_responses = []
    pulse_response_amps = []
    pulse_ids = []
    for pr in pair.pulse_responses:
        stim_pulse = pr.stim_pulse
        n_spikes = stim_pulse.n_spikes
        pulse_number = stim_pulse.pulse_number
        pulse_id = pr.stim_pulse_id
        ex_qc_pass = pr.ex_qc_pass
        in_qc_pass = pr.in_qc_pass
        pcr = stim_pulse.recording.patch_clamp_recording
        stim_freq = pcr.multi_patch_probe[0].induction_frequency
        clamp_mode = pcr.clamp_mode
        # current clamp
        if clamp_mode != 'ic':
            continue
        # ensure that there was only 1 presynaptic spike
        if n_spikes != 1:
            continue
        # we only want the first pulse of the train
        if pulse_number != 1:
            continue
        # only include frequencies up to 50Hz
        if stim_freq > 50:
            continue

        data = pr.data
        start_time = pr.start_time
        spike_time = stim_pulse.spikes[0].max_dvdt_time
        data_trace = Trace(data=data, t0=start_time - spike_time, sample_rate=db.default_sample_rate)

        if synapse_type == 'ex' and ex_qc_pass is True:
            pulse_responses.append(data_trace)
            pulse_ids.append(pulse_id)
            pulse_response_amps.append(pr.pulse_response_strength.pos_amp)
        if synapse_type == 'in' and in_qc_pass is True:
            pulse_responses.append(data_trace)
            pulse_ids.append(pulse_id)
            pulse_response_amps.append(pr.pulse_response_strength.neg_amp)

    return pulse_responses, pulse_ids, pulse_response_amps



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Populate first pulse table')
    parser.add_argument('--rebuild', action='store_true', default=False)
    parser.add_argument('--update', action='store_true', default=False)
    parser.add_argument('--limit', type=int)

    args = parser.parse_args(sys.argv[1:])

    if args.rebuild:
        args.rebuild = raw_input("Rebuild first pulse features? ") == 'y'
    if args.rebuild:
        first_pulse_features_tables.drop_tables()
        args.update = True


    init_tables()

    app = pg.mkQApp()
    pg.dbg()
    pg.setConfigOption('background', 'w')
    pg.setConfigOption('foreground', 'k')

    if args.update:
        update_analysis(limit=args.limit)





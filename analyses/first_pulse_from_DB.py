from multipatch_analysis.database import database as db
from neuroanalysis.data import Trace
from neuroanalysis.fitting import fit_psp
import argparse, sys, re
import pyqtgraph as pg

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
        'in_first_pulse_features': [
            """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed inhibitory qc in current clamp""",
            ('pair_id', 'pair.id', '', {'index': True}),
            ('in_ic_fit_amp', 'float', 'amplitude from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('in_ic_fit_latency', 'float', 'latency from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('in_ic_fit_rise_time', 'float', 'rise time from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('in_ic_fit_decay_tau', 'float', 'decay tau from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('in_ic_amp_cv', 'float', 'coefficient of variation for first pulse amplitude in 10, 20, 50 Hz stimuli'),
            ('in_avg_psp', 'array', 'average psp time series, spike aligned, baseline subtracted')
            #TODO: consider removing 50Hz responses from decay calculation
        ],

        'ex_first_pulse_features': [
            """Contains results of psp_fit on spike aligned, average first pulse PSP for each
            connection that passed excitatory qc in current clamp""",
            ('pair_id', 'pair.id', '', {'index': True}),
            ('ex_ic_fit_amp', 'float', 'amplitude from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ex_ic_fit_latency', 'float', 'latency from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ex_ic_fit_rise_time', 'float', 'rise time from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ex_ic_fit_decay_tau', 'float', 'decay tau from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ex_ic_amp_cv', 'float', 'coefficient of variation for first pulse amplitude in 10, 20, 50 Hz stimuli'),
            ('ex_avg_psp', 'array', 'average psp time series, spike aligned, baseline subtracted')
            # TODO: consider removing 50Hz responses from decay calculation
        ],

    }

    def create_mappings(self):
        TableGroup.create_mappings(self)

        INFirstPulseFeatures = self['in_first_pulse_features']
        EXFirstPulseFeatures = self['ex_first_pulse_features']

        db.Pair.in_first_pulse_features = db.relationship(INFirstPulseFeatures, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        INFirstPulseFeatures.pair = db.relationship(db.Pair, back_populates="in_first_pulse_features", single_parent=True)
        db.Pair.ex_first_pulse_features = db.relationship(EXFirstPulseFeatures, back_populates="pair", cascade="delete",
                                                          single_parent=True, uselist=False)
        EXFirstPulseFeatures.pair = db.relationship(db.Pair, back_populates="ex_first_pulse_features",
                                                    single_parent=True)


first_pulse_features_tables = FirstPulseFeaturesTableGroup()


def init_tables():
    global EXFirstPulseFeatures, INFirstPulseFeatures
    first_pulse_features_tables.create_tables()
    EXFirstPulseFeatures = first_pulse_features_tables['ex_first_pulse_features']
    INFirstPulseFeatures = first_pulse_features_tables['in_first_pulse_features']

def update_analysis(limit=None):
    s = db.Session()
    q = s.query(db.Pair, EXFirstPulseFeatures).outerjoin(EXFirstPulseFeatures).filter(EXFirstPulseFeatures.pair_id == None)
    if limit is not None:
        q = q.limit(limit)
    print("Updating %d pairs.." % q.count())
    records = q.all()
    for i,record in enumerate(records):
        pair = record[0]
        ex_pulse_responses, in_pulse_responses = filter_pulse_responses(pair)

        ex_results = first_pulse_features(ex_pulse_responses)
        in_results = first_pulse_features(in_pulse_responses)
        ex_fpf = EXFirstPulseFeatures(pair=pair, **ex_results) ## make keys exactly the column name
        in_fpf = INFirstPulseFeatures(pair=pair, **in_results)
        # fpf.pair = pair
        # fpf.ic_fit_amp = results['ic_fit_amp']
        # fpf.ic_fit_latency = results['fit_latency']
        # fpf.ic_fit_rise_time = results['fit_rise_time']
        # fpf.ic_fit_decay_tau = results['fit_decay_tau']
        # fpf.ic_amp_cv = results['amp_cv']
        # fpf.avg_psp = results['avg_psp_trace']
        # s.add(fpf)
        if i%10 == 0:
            s.commit()
    s.commit()
    s.close()


def first_pulse_features(pulse_responses):
    features = dict.fromkeys(['fit_amp', 'fit_latency', 'fit_rise_time', 'fit_decay_tau', 'amp_cv', 'avg_psp_trace'], 0)

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

    ex_pulse_responses = []
    in_pulse_responses = []
    for pr in pair.pulse_responses:
        stim_pulse = pr.stim_pulse
        n_spikes = stim_pulse.n_spikes
        pulse_number = stim_pulse.pulse_number
        ex_qc_pass = pr.ex_qc_pass
        in_qc_pass = pr.in_qc_pass
        pcr = stim_pulse.recording.patch_clamp_recording
        stim_freq = pcr.multi_patch_probe.induction_frequency
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
        if ex_qc_pass is True:
            ex_pulse_responses.append(data_trace)
        if in_qc_pass is True:
            in_pulse_responses.append(data_trace)

    pg.stack()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rebuild', action='store_true', default=False)
    parser.add_argument('--limit', type=int)
    args = parser.parse_args(sys.argv[1:])

    if args.rebuild:
        args.rebuild = raw_input("Rebuild first pulse features? ") == 'y'
    if args.rebuild:
        ex_first_pulse_features_tables.drop_tables()

    init_tables()

    update_analysis(limit=args.limit)

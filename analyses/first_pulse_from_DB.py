from multipatch_analysis.database import database as db
import argparse, sys

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
            connection in current clamp""",
            ('pair_id', 'pair.id', '', {'index': True}),
            ('ic_fit_amp', 'float', 'amplitude from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_latency', 'float', 'latency from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_rise_time', 'float', 'rise time from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_decay_tau', 'float', 'decay tau from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_amp_cv', 'float', 'coefficient of variation for first pulse amplitude in 10, 20, 50 Hz stimuli'),
            ('avg_psp', 'array', 'average psp time series, spike aligned, baseline subtracted')
            #TODO: consider removing 50Hz responses from decay calculation
        ],
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
        results = first_pulse_features(pair)
        fpf = FirstPulseFeatures()
        fpf.pair = pair
        fpf.ic_fit_amp = results['fit_amp']
        s.add(fpf)
        if i%10 == 0:
            s.commit()
    s.commit()
    s.close()


def first_pulse_features(pair):
    return {'fit_amp': 0}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rebuild', action='store_true', default=False)
    parser.add_argument('--limit', type=int)
    args = parser.parse_args(sys.argv[1:])

    if args.rebuild:
        args.rebuild = raw_input("Rebuild first pulse features? ") == 'y'
    if args.rebuild:
        first_pulse_features_tables.drop_tables()

    init_tables()

    update_analysis(limit=args.limit)

from multipatch_analysis.database import database as db

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
            ('pair_id', 'pair.id', '', {'index': True}),
            ('ic_fit_amp', 'float', 'amplitude from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_latency', 'float', 'latency from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_rise_time', 'float', 'rise time from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_fit_decay_tau', 'float', 'decay tau from psp_fit to first pulse avg of 10, 20, 50 Hz stimuli'),
            ('ic_amp_cv', 'float', 'coefficient of variatiation for first pulse amplitude in 10, 20, 50 Hz stimuli'),
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


if __name__ == '__main__':
    init_tables()
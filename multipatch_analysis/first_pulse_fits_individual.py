"""seperating individal fitting from average and instead querying the database for the average fits.
"""

#import pdb
from neuroanalysis.data import Trace, TraceList
from multipatch_analysis.database import database as db
import multipatch_analysis.connection_strength as cs 
from multipatch_analysis.database.database import TableGroup
import matplotlib.pyplot as plt
import numpy as np
import time
from neuroanalysis.fitting import fit_psp

class IndFirstPulseFitTableGroup(TableGroup):
    """Fits first pulse for each individual sweeps.
    """
    schemas = {
        'individual_first_pulse_fit': [
            """Best parameters fit to individual first pulses with initial conditions set 
            by the fit of the average first pulse available in the avg_first_pulse_fit_force_sign table.""",
            ('pulse_response_id', 'pulse_response.id', '', {'index': True, 'unique': True})].append(common)
    }
    
    def create_mappings(self):
        TableGroup.create_mappings(self)
        IndividualFirstPulseFits = self['individual_first_pulse_fit']

        here might want to back populate by pulse response not pair
        db.Pair.individual_first_pulse_fit = db.relationship(IndividualFirstPulseFits, back_populates="pair", cascade="delete",
                                                      single_parent=True, uselist=False)
        IndividualFirstPulseFits.pair = db.relationship(db.Pair, back_populates="individual_first_pulse_fit", single_parent=True)

ind_first_pulse_fit_table = IndFirstPulseFitTableGroup()

def init_tables():
    global IndividualFirstPulseFits 
    ind_first_pulse_fit_table.create_tables()

    IndividualFirstPulseFits = ind_first_pulse_fit_table['individual_first_pulse_fit']

# create tables in database and add global variables for ORM classes
init_tables()
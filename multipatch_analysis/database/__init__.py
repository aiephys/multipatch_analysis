from .database import Session, aliased, default_session, get_default_session, reset_db, vacuum, dispose_engines, default_sample_rate, db_name, bake_sqlite

# Import table definitions from DB modules
from .pipeline import *
from .slice import *
from .experiment import *
from .morphology import *
from .dataset import *
from .pulse_response_strength import *
from .dynamics import *
from .connection_strength import *
from .avg_first_pulse_fit import *
from .single_first_pulse_fit import *
from .gap_junction import *


@default_session
def slice_from_timestamp(ts, session=None):
    slices = session.query(Slice).filter(Slice.acq_timestamp==ts).all()
    if len(slices) == 0:
        raise KeyError("No slice found for timestamp %0.3f" % ts)
    elif len(slices) > 1:
        raise KeyError("Multiple slices found for timestamp %0.3f" % ts)
    
    return slices[0]


@default_session
def experiment_from_timestamp(ts, session=None):
    expts = session.query(Experiment).filter(Experiment.acq_timestamp==ts).all()
    if len(expts) == 0:
        # For backward compatibility, check for timestamp truncated to 2 decimal places
        for expt in session.query(Experiment).all():
            if abs((expt.acq_timestamp - ts)) < 0.01:
                return expt
        
        raise KeyError("No experiment found for timestamp %0.3f" % ts)
    elif len(expts) > 1:
        raise RuntimeError("Multiple experiments found for timestamp %0.3f" % ts)
    
    return expts[0]


@default_session
def list_experiments(session=None):
    return session.query(Experiment).all()


def query(*args, **kwds):
    return get_default_session().query(*args, **kwds)


def rollback():
    get_default_session().rollback()



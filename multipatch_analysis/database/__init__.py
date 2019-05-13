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


def query_pairs(project_name=None, acsf=None, age=None, species=None, distance=None, session=None, internal=None):
    """Generate a query for selecting pairs from the database.

    Parameters
    ----------
    project_name : str
        Value to match from experiment.project_name (e.g. "mouse V1 coarse matrix" or "human coarse matrix")
    """
    pre_cell = aliased(Cell, name='pre_cell')
    post_cell = aliased(Cell, name='post_cell')
    pairs = session.query(
        Pair, 
        # pre_cell,
        # post_cell,
        # Experiment,
        # Pair.synapse,
        # ConnectionStrength.synapse_type,
    )
    pairs = pairs.join(pre_cell, pre_cell.id==Pair.pre_cell_id).join(post_cell, post_cell.id==Pair.post_cell_id)
    pairs = pairs.join(Experiment).join(Slice)
    pairs = pairs.join(ConnectionStrength)
    
    if project_name is not None:
        if isinstance(project_name, str):
            pairs = pairs.filter(Experiment.project_name==project_name)
        else:
            pairs = pairs.filter(Experiment.project_name.in_(project_name))

    if acsf is not None:
        if isinstance(acsf, str):
            pairs = pairs.filter(Experiment.acsf==acsf)
        else:
            pairs = pairs.filter(Experiment.acsf.in_(acsf))

    if age is not None:
        if age[0] is not None:
            pairs = pairs.filter(Slice.age>=age[0])
        if age[1] is not None:
            pairs = pairs.filter(Slice.age<=age[1])

    if distance is not None:
        if distance[0] is not None:
            pairs = pairs.filter(Pair.distance>=distance[0])
        if distance[1] is not None:
            pairs = pairs.filter(Pair.distance<=distance[1])

    if species is not None:
        pairs = pairs.filter(Slice.species==species)

    if internal is not None:
        if isinstance(internal, str):
            pairs = pairs.filter(Experiment.internal==internal)
        else:
            pairs = pairs.filter(Experiment.internal.in_(internal))

    return pairs

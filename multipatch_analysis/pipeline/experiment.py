import os, glob, re, time
import numpy as np
from datetime import datetime
from collections import OrderedDict
from acq4.util.DataManager import getDirHandle
from .. import database as db
from ..database import experiment_tables
from .pipeline_module import DatabasePipelineModule
from .. import config, synphys_cache
from .. import lims
from ..util import datetime_to_timestamp
from ..experiment import Experiment
from .slice import SlicePipelineModule


class ExperimentPipelineModule(DatabasePipelineModule):
    """Imports per-experiment metadata into DB.
    """
    name = 'experiment'
    dependencies = [SlicePipelineModule]
    table_group = experiment_tables
    
    
    @classmethod
    def process_job(cls, job_id):
        # slices = all_slices()
        # path = slices[job_id]
        # dh = getDirHandle(path)
        # info = dh.info()
        

        # sl = db.Slice(**fields)
        # # sl.meta = {'db_timestamp': time.time()}
        cache = synphys_cache.get_cache()
        all_expts = cache.list_experiments()
        site_path = all_expts[job_id]
        expt = Experiment(site_path=site_path)
        
        session = db.Session(readonly=False)
        # try:
        #     session.add(sl)
        #     session.commit()
        # finally:
        #     session.close()
        
        # look up slice record in DB
        ts = expt.slice_timestamp
        slice_entry = db.slice_from_timestamp(ts, session=session)
        
        expt_info = expt.expt_info
        fields = {
            'original_path': expt.original_path,
            'storage_path': expt.server_path,
            'ephys_file': None if expt.nwb_file is None else os.path.relpath(expt.nwb_file, expt.path),
            'rig_name': expt.rig_name,
            'project_name': expt.project_name,
            'acq_timestamp': expt.timestamp,
            'target_region': expt_info.get('region'),
            'internal': expt_info.get('internal'),
            'acsf': expt_info.get('solution'),
            'target_temperature': expt.target_temperature,
        }

        # Create entry in experiment table
        expt_entry = db.Experiment(**fields)
        expt_entry.slice = slice_entry
        session.add(expt_entry)

        # create pipette and cell entries
        cell_entries = {}
        for e_id, elec in expt.electrodes.items():
            elec_entry = db.Electrode(experiment=expt_entry, ext_id=elec.electrode_id, device_id=elec.device_id)
            for k in ['patch_status', 'start_time', 'stop_time',  
                      'initial_resistance', 'initial_current', 'pipette_offset',
                      'final_resistance', 'final_current']:
                if hasattr(elec, k):
                    setattr(elec_entry, k, getattr(elec, k))
            session.add(elec_entry)

            if elec.cell is not None:
                cell = elec.cell
                cell_entry = db.Cell(
                    electrode=elec_entry,
                    ext_id=cell.cell_id,
                    cre_type=cell.cre_type,
                    target_layer=cell.target_layer,
                    is_excitatory=cell.is_excitatory,
                    depth=cell.depth,
                    position=cell.position,
                )
                session.add(cell_entry)
                cell_entries[cell] = cell_entry

        # create pairs
        for i, pre_cell in expt.cells.items():
            for j, post_cell in expt.cells.items():
                if i == j:
                    continue
                # check to see if i,j is in manual connection calls
                # (do not use expt.connections, which excludes some connections based on QC)
                syn_calls = expt.connection_calls
                synapse = None if syn_calls is None else ((i, j) in syn_calls)
                gap_calls = expt.gap_calls
                electrical = None if gap_calls is None else ((i, j) in gap_calls)

                pre_cell_entry = cell_entries[pre_cell]
                post_cell_entry = cell_entries[post_cell]
                p1, p2 = pre_cell.position, post_cell.position
                if None in [p1, p2]:
                    distance = None
                else:
                    distance = np.linalg.norm(np.array(p1) - np.array(p2))
                
                pair_entry = db.Pair(
                    experiment=expt_entry,
                    pre_cell=pre_cell_entry,
                    post_cell=post_cell_entry,
                    synapse=synapse,
                    electrical=electrical,
                    n_ex_test_spikes=0,  # will be counted later
                    n_in_test_spikes=0,
                    distance=distance,
                )
                session.add(pair_entry)

                pre_id = pre_cell_entry.electrode.device_id
                post_id = post_cell_entry.electrode.device_id

        expt_entry.meta = {'experiment_finished': time.time()}

        try:
            session.commit()
        except:
            session.rollback()
            raise
        finally:
            session.close()

    @classmethod
    def drop_jobs(cls, job_ids):
        """Remove all results previously stored for a list of job IDs.
        """
        session = db.Session(readonly=False)
        expts = session.query(db.Experiment).filter(db.Experiment.acq_timestamp.in_(job_ids))
        for exp in expts:
            session.delete(exp)
        session.commit()
        
    @classmethod
    def finished_jobs(cls):
        """Return an ordered dict of job IDs that have been processed by this module and
        the dates when they were processed.

        Note that some results returned may be obsolete if dependencies have changed.
        """
        session = db.Session()
        expts = session.query(db.Experiment.acq_timestamp, db.Experiment.time_created).all()
        session.rollback()
        return OrderedDict([(uid, datetime_to_timestamp(date)) for uid, date in expts])

    @classmethod
    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        finished_slices = SlicePipelineModule.finished_jobs()
        
        cache = synphys_cache.get_cache()
        all_expts = cache.list_experiments()
        
        n_errors = 0
        ready = OrderedDict()
        for ts, path in all_expts.items():
            try:
                expt = Experiment(site_path=path)
            except Exception:
                n_errors += 1
                continue
            raw_data_mtime = expt.last_modification_time
            slice_ts = expt.slice_timestamp
            slice_mtime = finished_slices.get(slice_ts, None)
            if slice_mtime is None:
                continue
            ready[ts] = max(raw_data_mtime, slice_mtime)
        
        print("Found %d experiments; %d are able to be processed, %d were skipped due to errors." % (len(all_expts), len(ready), n_errors))
        return ready

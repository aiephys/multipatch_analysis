from aisynphys.pipeline.pipeline_module import DatabasePipelineModule
#import multipatch_analysis.database as db
from aisynphys import config
from .opto_slice import OptoSlicePipelineModule
from collections import OrderedDict
import csv, codecs, glob, os
from acq4.util.DataManager import getDirHandle
from neuroanalysis.data.experiment import Experiment
from neuroanalysis.data.libraries import opto
from . import data_model

##### TODO: GO BACK TO EXPERIEMENT BEING DEPENDENT ON SLICE -- IN ALL SLICES USE EXPERIMENTS.CSV TO COME UP WITH SLICE LIST

class OptoExperimentPipelineModule(DatabasePipelineModule):
    """Imports per-experiment metadata into DB.
    """
    name = 'opto_experiment'
    dependencies = [OptoSlicePipelineModule]
    #table_group = db.experiment_tables
    table_group = ['experiment', 'electrode', 'cell', 'pair']

    @classmethod
    def create_db_entries(cls, job, session, expt=None):
        """Generate DB entries for *job_id* and add them to *session*.
        """
        job_id = job['job_id']
        db = job['database']

        try:
            if expt is None:
                expt = load_experiment(job_id)
            # #    raise Exception("Please pass an Experiment (data model) object. Looking up job_ids is not yet implemented.")
            #     all_expts = read_expt_csvs()
            #     indices = [i for i, e in enumerate(all_expts['expt_list']) if job_id in e['experiment']]
            #     if len(indices) > 1:
            #         raise Exception("Cannot resolve job_id: %s. Found %s" % (job_id, [all_expts['expt_list'][i]['experiment'] for i in indices]))
            #     elif len(indices) == 0:
            #         raise Exception("Could not find csv entry for %s"%job_id)

            #     entry = all_expts['expt_list'][indices[0]]
            #     entry['distances'] = [e for e in all_expts['distances'] if e['exp_id']==job_id]
            #     #print('create_db_entries for:', entry['site_path'], "job_id:", job_id)
            #     if entry['site_path'] != '':
            #         expt = Experiment(site_path=entry['site_path'], loading_library=opto, meta_info=entry)
            #     else:
            #         cnx_json = os.path.join(config.connections_dir, entry['experiment'])
            #         expt = Experiment(load_file=cnx_json, loading_library=opto, meta_info=entry)

            # look up slice record in DB
            try:
                ts = expt.slice_timestamp
            except KeyError:
                ts = 0.0
                

            slice_entry = db.slice_from_timestamp(ts, session=session)

            expt_info = expt.expt_info
            if expt_info is None:
                expt_info = {}

            # dig to find out which rig this was recorded on
            rig = expt_info.get('rig_name', None)
            if rig is None and expt.data is not None:
                ## serial number is recorded in many places, make sure they converge on one rig
                sns = []
                for sweeps in expt.data.notebook().values():
                    for channel in sweeps:
                        sn = channel.get('Serial Number', None)
                        if sn is not None:
                            sns.append(sn)
                unique_sns = list(set(sns))
                rigs = []
                for sn in unique_sns:
                    rigs.append(data_model.get_rig_name_from_serial_number(sn))
                unique_rigs = list(set(rigs))
                if len(unique_rigs) != 1:
                    raise Exception("Could not resolve rig for experiment %s. Found %s" %(expt.uid, unique_rigs))
                rig = unique_rigs[0]

                #serial_number = expt.nwb.notebook()[0][0]['Serial Number']
                #rig = data_model.get_rig_name_from_serial_number(serial_number)



            fields = {
                'storage_path': expt.original_path, 
                'ephys_file': None if expt.ephys_file is None else os.path.relpath(expt.ephys_file, expt.path),
                'rig_name': rig,
                #'project_name': expt.project_name,
                'acq_timestamp': expt.timestamp,
                #'target_region': expt_info.get('region'),
                'internal': expt_info.get('internal'),
                'acsf': expt_info.get('solution'),
                #'target_temperature': expt.target_temperature,
                'ext_id': expt.uid
            }

            expt_entry = db.Experiment(**fields)
            expt_entry.slice = slice_entry
            session.add(expt_entry)

            cell_entries = {}
            for name, cell in expt.cells.items():
                if cell.electrode is not None:
                    elec = cell.electrode
                    elec_entry = db.Electrode(experiment=expt_entry, ext_id=elec.electrode_id, device_id=elec.device_id)
                    for k in ['patch_status', 'start_time', 'stop_time',  
                        'initial_resistance', 'initial_current', 'pipette_offset',
                        'final_resistance', 'final_current']:
                        if hasattr(elec, k):
                            setattr(elec_entry, k, getattr(elec, k))
                    session.add(elec_entry)

                cell_entry = db.Cell(
                    experiment=expt_entry,
                    electrode=elec_entry if cell.electrode is not None else None,
                    ext_id=cell.cell_id,
                    cre_type=cell.cre_type,
                    target_layer=cell.target_layer,
                    is_excitatory=cell.is_excitatory,
                    depth=cell.depth,
                    position=cell.position,
                )
                session.add(cell_entry)
                cell_entries[cell] = cell_entry

            for name, pair in expt.pairs.items():
                if pair._connection_call == 'excitatory':
                    sign = +1
                elif pair._connection_call == 'inhibitory':
                    sign = -1
                else:
                    sign = 0
                pair_entry = db.Pair(
                    experiment=expt_entry,
                    pre_cell=cell_entries[pair.preCell],
                    post_cell=cell_entries[pair.postCell],
                    has_synapse=pair.isSynapse(),
                    has_electrical=None,
                    n_ex_test_spikes=0,  # will be counted in opto_dataset_pipeline_module
                    n_in_test_spikes=0,
                    distance=pair.distance,
                    #synapse_sign = sign
                )
                session.add(pair_entry)

            #session.commit()

        except:
            session.rollback()
            raise
        
    
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        # only need to return from experiment table; other tables will be dropped automatically.
        db = self.database
        return session.query(db.Experiment).filter(db.Experiment.ext_id.in_(job_ids)).all()

    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """

        slice_module = self.pipeline.get_module('opto_slice')
        finished_slices = slice_module.finished_jobs()
        
        # cache = synphys_cache.get_cache()
        # all_expts = cache.list_experiments()
        db = self.database
        session = db.session()
        slices = session.query(db.Slice.storage_path).all()
        slice_paths = [s[0] for s in slices]
        
        #ymls = []
        #for rec in slices:
        #    path = rec[0]
        #    ymls.extend(glob.glob(os.path.join(config.synphys_data, path, 'site_*', 'pipettes.yml')))
        expts = read_expt_csvs()
        
        n_errors = 0
        n_no_slice = 0
        ready = OrderedDict()
        print("checking for ready expts....")
        for i, expt in enumerate(expts['expt_list']):
            #print("Checking experiment %i/%i"%(i, len(expts['expt_list'])))
            site_path = os.path.join(config.synphys_data, expt['site_path'])
            slice_path = getDirHandle(os.path.split(site_path)[0]).name(relativeTo=getDirHandle(config.synphys_data))
            #print slice_paths
            #if not slice_path in slice_paths:
            #    #print("Did not find slice path for %s"%slice_path)
            #    n_no_slice += 1
            #    continue
            try:
                if expt['site_path'] == '':
                    cnx_json = os.path.join(config.connections_dir, expt['experiment'])
                    ex = Experiment(load_file=cnx_json, loading_library=opto, meta_info=expt)
                else:
                    ex = Experiment(site_path=site_path, loading_library=opto)

                raw_data_mtime = ex.last_modification_time
                try:
                    slice_ts = ex.slice_timestamp
                except KeyError:
                    slice_ts = 0.0
                slice_mtime, slice_success = finished_slices.get('%.3f'%slice_ts, (None, None))
                #print('found expt for path:', site_path)
            except Exception:
                n_errors += 1
                continue
            if slice_mtime is None or slice_success is False:
            #    slice_mtime = 0
                n_no_slice += 1
                continue

            ready[ex.uid] = max(raw_data_mtime, slice_mtime)
        
        print("Found %d experiments; %d are able to be processed, %d were skipped due to errors, %d were skipped due to missing or failed slice entries." % (len(expts['expt_list']), len(ready), n_errors, n_no_slice))
        return ready

    def dependent_job_ids(self, module, job_ids):
        """Return a list of all finished job IDs in this module that depend on 
        specific jobs from another module.
        """
        if type(module) not in self.dependencies:
            raise ValueError("%s does not depend on module %s" % (cls, module))
        
        db = self.database
        session = db.session()
        dep_ids = session.query(db.Experiment.ext_id).join(db.Slice).filter(db.Slice.acq_timestamp.in_(job_ids)).all()
        session.rollback()
        return [rec.ext_id for rec in dep_ids]


_csv_data = None
def read_expt_csvs():

    global _csv_data
    if _csv_data is not None:
        return _csv_data

    _csv_data = OrderedDict()

    expt_csv = config.experiment_csv
    distance_csv = config.distance_csv

    _csv_data['expt_list'] = []
    with open(expt_csv, 'r') as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            _csv_data['expt_list'].append(row)

    _csv_data['distances']=[]
    if distance_csv is not None:
        with codecs.open(distance_csv, "r", encoding="utf-8-sig") as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader:
                _csv_data['distances'].append(row)

    return _csv_data

def load_experiment(job_id):
    all_expts = read_expt_csvs()
    indices = [i for i, e in enumerate(all_expts['expt_list']) if job_id in e['experiment']]
    if len(indices) > 1:
        raise Exception("Cannot resolve job_id: %s. Found %s" % (job_id, [all_expts['expt_list'][i]['experiment'] for i in indices]))
    elif len(indices) == 0:
        raise Exception("Could not find csv entry for %s"%job_id)

    entry = all_expts['expt_list'][indices[0]]
    entry['distances'] = [e for e in all_expts['distances'] if e['exp_id']==job_id]
    #print('create_db_entries for:', entry['site_path'], "job_id:", job_id)
    if entry['site_path'] != '':
        expt = Experiment(site_path=os.path.join(config.synphys_data, entry['site_path']), loading_library=opto, meta_info=entry)
    else:
        cnx_json = os.path.join(config.connections_dir, entry['experiment'])
        expt = Experiment(load_file=cnx_json, loading_library=opto, meta_info=entry)

    return expt


# coding: utf8
"""
For generating a DB table describing cell morphology.

"""
from __future__ import print_function, division

import os, datetime, hashlib, json
from sqlalchemy.orm import joinedload
from collections import OrderedDict
from ...util import timestamp_to_datetime, optional_import
from ...data.pipette_metadata import PipetteMetadata
from ... import config, lims
from .pipeline_module import MultipatchPipelineModule
from .experiment import ExperimentPipelineModule


col_names = {
    'Qual_Morpho_Type': {'name': 'qual_morpho_type', 'type': 'str'},
    'dendrite_type': {'name': 'dendrite_type', 'type': ['spiny', 'aspiny', 'sparsely spiny', 'NEI']},
    'Apical_truncation_distance': {'name': 'apical_trunc_distance', 'type': 'float', 'scale': 1e-6},
    'Axon_truncation_distance': {'name': 'axon_trunc_distance', 'type': 'float', 'scale': 1e-6},
    'Axon origination': {'name': 'axon_origin', 'type': ['soma', 'dendrite', 'unclear', 'NEI']},
    'Axon_truncation': {'name': 'axon_truncation', 'type': ['truncated', 'borderline', 'intact', 'unclear', 'NEI']},
    'Apical_truncation': {'name': 'apical_truncation', 'type': ['truncated', 'borderline', 'intact','unclear', 'NEI']},
}


class MorphologyPipelineModule(MultipatchPipelineModule):
    """Imports cell morphology data for each experiment
    """
    name = 'morphology'
    dependencies = [ExperimentPipelineModule]
    table_group = ['morphology']
    
    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        job_id = job['job_id']

        # Load experiment from DB
        expt = db.experiment_from_timestamp(job_id, session=session)
        morpho_results = morpho_db()
        
        path = os.path.join(config.synphys_data, expt.storage_path)
        pip_meta = PipetteMetadata(path)

        for cell_id,cell in expt.cells.items():
            # How the experimenter described the morphology
            user_morpho = pip_meta.pipettes[cell.ext_id].get('morphology')
            cell_specimen_id = cell.meta.get('lims_specimen_id')
            cell_morpho = morpho_results.get(cell_specimen_id)

            if user_morpho in (None, ''):
                pyramidal = None
            elif user_morpho == 'pyr':
                pyramidal = True
            else:
                print("Unknown morphology string: %s" % user_morpho)
                pyramidal = None

            results = {
                'pyramidal': pyramidal,
            }
            
            morpho_db_hash = hash_record([cell_morpho])
            results['meta'] = {'morpho_db_hash': morpho_db_hash}
            
            if cell_morpho is not None:  
                for morpho_db_name, result in col_names.items():
                    col_name = result['name']
                    col_type = result['type']
                    data = cell_morpho[morpho_db_name]
                    if data is not None:
                        if isinstance(col_type, list):
                            d = [t for t in col_type if t == data]
                            if len(d) == 1:
                                data = d[0]
                            else:
                                raise Exception ('Error parsing morphology annotation %s for cell %d from column %s' % (data, cell_specimen_id, morpho_db_name))
                        elif col_type == 'float':
                            try:
                                scale = result.get('scale', 1)
                                data = float(data) * scale
                            except ValueError:
                                raise Exception ('Error parsing morphology annotation %s for cell %d from column %s' % (data, cell_specimen_id, morpho_db_name))
                    results[col_name] = data

            # Write new record to DB
            morphology = db.Morphology(cell_id=cell.id, **results)
            
            # Update cell_class_nonsynaptic
            #  (cell_class gets updated later)
            dtype = results.get('dendrite_type', None)
            morpho_class = {'spiny': 'ex', 'aspiny': 'in', 'sparsely spiny': 'in'}.get(dtype, None)
            cell_meta = cell.meta.copy()
            cell_meta['morpho_cell_class'] = morpho_class
            cell.meta = cell_meta

            # this gets updated again in later modules
            cell.cell_class, cell.cell_class_nonsynaptic = cell._infer_cell_classes()
                
            session.add(morphology)
        
    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        return session.query(db.Morphology).filter(db.Morphology.cell_id==db.Cell.id).filter(db.Cell.experiment_id==db.Experiment.id).filter(db.Experiment.ext_id.in_(job_ids)).all()

    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        db = self.database
        # All experiments and their creation times in the DB
        expts = self.pipeline.get_module('experiment').finished_jobs()

        # Look up nwb file locations for all experiments
        session = db.session()
        # expts = session.query(db.Experiment).filter(db.Experiment.ext_id==1521667891.153).all()
        session.rollback()
        
        # Return the greater of NWB mod time and experiment DB record mtime
        ready = OrderedDict()

        try:
            morpho_results = morpho_db()
        except ImportError as exc:
            print("Skipping morphology: %s" % str(exc))
            return ready
            
        for expt_id, (expt_mtime, success) in expts.items():
            if success is not True:
                continue

            q = session.query(db.Experiment)
            # joinedload should speed up access to cell/morpho attributes by eager loading along with the experiment
            q = q.options(joinedload(db.Experiment.cell_list).joinedload(db.Cell.morphology))
            q = q.filter(db.Experiment.ext_id==expt_id)
            expt = q.all()[0]

            ready[expt_id] = {'dep_time': expt_mtime}
            needs_update = False
            for cell_ext_id, cell in expt.cells.items():
                if cell.morphology is None:
                    needs_update = True
                    break
                cell_specimen_id = cell.meta.get('lims_specimen_id')
                morpho_rec = morpho_results.get(cell_specimen_id, None)
                morpho_db_hash = hash_record([morpho_rec])
                prev_hash = None if cell.morphology.meta is None else cell.morphology.meta['morpho_db_hash']
                if morpho_db_hash != prev_hash:
                    needs_update = True
                    break

            if needs_update:        
                ready[expt_id] = {'dep_time': datetime.datetime.now()}
        
        return ready


def hash_record(rec):
    return hashlib.md5(repr(rec).encode()).hexdigest()

def import_morpho_db():
    import pyodbc
    cnxn_str = r'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)}; DBQ=%s' % config.morpho_address
    cnxn = pyodbc.connect(cnxn_str)
    cursor = cnxn.cursor()
    morpho_table = cursor.execute('select * from MPATCH_CellsofCluster')
    results = morpho_table.fetchall()
    fields = [r[0] for r in results[0].cursor_description]
    morpho_results = {int(r.cell_specimen_id): {k:getattr(r, k) for k in fields} for r in results}
    
    return morpho_results

morpho_cache = None
def morpho_db():
    global morpho_cache
    if morpho_cache is None:
        if hasattr(config, 'morpho_address'):
            morpho_cache = import_morpho_db()
        else:
            # json requires string keys, so we have to convert back to int here:
            morpho_cache = {int(k):v for k,v in json.load(open(config.morpho_json_file)).items()}
    
    return morpho_cache

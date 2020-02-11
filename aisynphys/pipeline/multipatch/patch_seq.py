# coding: utf8
"""
For generating a DB table describing cell transcriptomics.

"""
from __future__ import print_function, division

import os, datetime, shutil, tempfile, hashlib
import pandas as pd
from collections import OrderedDict
from ...util import timestamp_to_datetime, optional_import
from acq4.util.DataManager import getHandle
from ... import config
from .pipeline_module import MultipatchPipelineModule
from .experiment import ExperimentPipelineModule
pyodbc = optional_import('pyodbc')

amp_cols = {
            'Comment': 'meta',
            'Result pass/fail BA': 'result_BA',
            '% area 400-10000bp BA': 'area_400_10000bp',
            'Picogreen pg/ul': 'picogreen_yield',
            }    

mapping_cols = {
            'cluster_detail': 'cluster_detail',
            'cluster_label':'cluster_label',
            'score': 'score',
            'res_index': 'res_index',
            'topLeaf': 'top_leaf',
            'topLeafValue': 'top_leaf_score',
            'broad_class_label': 'broad_class_label',
            'subclass_label': 'sublass_label',
            'quality_score_label': 'quality_score',
            'seurat_cluster_label': 'seurat_cluster',
            'seurat_prediction_score_label': 'seurat_score',
            'Tree_first_cl': 'tree_first_cluster',
            'Tree_first_bt': 'tree_first_score',
            'Tree_second_cl': 'tree_second_cluster',
            'Tree_second_bt': 'tree_second_score',
            'Tree_third_cl': 'tree_third_cluster',
            'Tree_third_bt': 'tree_third_score',
            'Tree_call': 'tree_call',
            'Genes.Detected.CPM': 'genes_detected',
            'marker_sum_norm_label': 'norm_marker_sum'
            }

col_names = amp_cols.copy()
col_names.update(mapping_cols)

nucleus_bool = {'+': True, '-': False, '': None}

class PatchSeqPipelineModule(MultipatchPipelineModule):
    """Imports transcriptomic data for cells in patchseq experiments
    """
    name = 'patch_seq'
    dependencies = [ExperimentPipelineModule]
    table_group = ['patch_seq']

    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        job_id = job['job_id']

        # Load experiment from DB
        expt = db.experiment_from_timestamp(job_id, session=session)
        amp_results = get_amp_results()
        mapping_results = get_mapping_results()

        path = os.path.join(config.synphys_data, expt.storage_path)
        site_info = getHandle(path).info()
        headstages = site_info.get('headstages')
        if headstages is not None:
            patchseq_tubes = {hs_name.split('HS')[1]: hs['Tube ID'] for hs_name, hs in headstages.items()}
            nucleus = {hs_name.split('HS')[1]: hs['Nucleus'] for hs_name, hs in headstages.items()}
            no_tubes = all(t == '' for t in patchseq_tubes.values())
            if no_tubes is False:

                for cell_ext_id, cell in expt.cells.items():
                    tube_id = patchseq_tubes.get(cell_ext_id, '').strip()
                    if tube_id == '':
                        continue

                    results = {
                        'tube_id': tube_id,
                        'nucleus': nucleus_bool[nucleus.get(cell_ext_id, '')],
                        'patchseq_hash': None,
                    }

                    amp_result = amp_results.get(tube_id, {})
                    mapping_result = mapping_results.get(tube_id, {})
    
                    patchseq_results = amp_result.copy()
                    patchseq_results.update(mapping_result)
                    if patchseq_results is None:
                        continue

                    patchseq_hash = hashlib.md5(str(tuple(patchseq_results.values())).encode()).hexdigest()
                    results['patchseq_hash'] = patchseq_hash

                    for result_name, col_name in col_names.items():
                        data = patchseq_results.get(result_name)
                        if data is not None:
                            if col_name == 'meta':
                                data = {'amplification_comments': data}
                            if col_name == 'genes_detected':
                                data = int(data)
                            results[col_name] = data

                    tree_call = results.get('tree_call')
                    if tree_call is not None and tree_call in ['Core', 'I1']:
                        results['t_type'] = results['tree_first_cluster']

                    # Write new record to DB
                    patch_seq = db.PatchSeq(cell_id=cell.id, **results)
                    session.add(patch_seq)

    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
            
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        db = self.database
        return session.query(db.PatchSeq).filter(db.PatchSeq.cell_id==db.Cell.id).filter(db.Cell.experiment_id==db.Experiment.id).filter(db.Experiment.ext_id.in_(job_ids)).all()

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
            amp_results = get_amp_results()
            mapping_results = get_mapping_results()
        except ImportError as exc:
            print("Skipping patchseq: %s" % str(exc))
            return ready

        patchseq_results = amp_results.copy()
        patchseq_results.update(mapping_results)

        for expt_id, (expt_mtime, success) in expts.items():
            if success is not True:
                continue

            expt = session.query(db.Experiment).filter(db.Experiment.ext_id==expt_id).all()[0]
            ready[expt_id] = {'dep_time': expt_mtime}

            path = os.path.join(config.synphys_data, expt.storage_path)
            site_info = getHandle(path).info()
            headstages = site_info.get('headstages')
            if headstages is None:
                continue

            patchseq_tubes = {hs_name.split('HS')[1]: hs['Tube ID'] for hs_name, hs in headstages.items()}
            if patchseq_tubes is None:
                continue

            patchseq_hash_compare = []
            for cell_ext_id, cell in expt.cells.items():
                tube_id = patchseq_tubes.get(cell_ext_id, '').strip()
                if tube_id not in patchseq_results:
                    continue
                patchseq_data = patchseq_results[tube_id]
                patchseq_data_hash = hashlib.md5(str(tuple(patchseq_data.values())).encode()).hexdigest()
                patchseq_rec = session.query(db.PatchSeq).join(db.Cell).filter(db.Cell.id == cell.id).all()
                if len(patchseq_rec) == 1:
                    patchseq_rec_hash = patchseq_rec[0].patchseq_hash
                    patchseq_hash_compare.append(patchseq_data_hash == patchseq_rec_hash)
                else:
                    patchseq_hash_compare.append(False)
            if all(patchseq_hash_compare) is False:
                ready[expt_id] = {'dep_time': datetime.datetime.now()}

        return ready

amp_cache = None
def get_amp_results(): 
    global amp_cache
    if amp_cache is None:
        folder_path = config.amplification_report_address
        file_paths = [os.path.join(folder_path, path) for path in os.listdir(folder_path)]
        mod_time = {os.path.getmtime(path): path for path in file_paths}
        most_recent = max(mod_time.keys())
        current_results = mod_time[most_recent]
        open_file, local_copy = tempfile.mkstemp(prefix='amplification_results', suffix='.xlsx')
        shutil.copy(current_results, local_copy)
        amp_results = pd.read_excel(local_copy, header=2, index_col=False)
        amp_results.set_index('Sample ID', inplace=True)
        amp_results_reduced = amp_results[list(amp_cols.keys())]
        amp_cache = amp_results_reduced.to_dict('index')

        os.close(open_file)
        os.remove(local_copy)

    return amp_cache

mapping_cache = None
def get_mapping_results():
    global mapping_cache
    if mapping_cache is None:
        folder_path = config.mapping_report_address
        mouse_path = folder_path + '/mouse_patchseq_VISp_current/mapping.df.with.bp.40.lastmap.csv'
        human_path = folder_path + '/human/human_patchseq_MTG_current/mapping.df.lastmap.csv'
        mouse_open_file, mouse_local = tempfile.mkstemp(prefix='mapping_mouse', suffix='.csv')
        shutil.copy(mouse_path, mouse_local)
        mouse_results = pd.read_csv(mouse_local, header=0, index_col=False, dtype=object)
        human_open_file, human_local = tempfile.mkstemp(prefix='mapping_human', suffix='.csv')
        shutil.copy(human_path, human_local)
        human_results = pd.read_csv(human_local, header=0, index_col=False, dtype=object)

        mapping_results = pd.concat([mouse_results, human_results], ignore_index=True, sort=False)
        mapping_results.set_index('sample_id', inplace=True)
        mapping_results_reduced = mapping_results[list(mapping_cols.keys())]
        mapping_cache = mapping_results.to_dict('index')

        os.close(mouse_open_file)
        os.remove(mouse_local)
        os.close(human_open_file)
        os.remove(human_local)

    return mapping_cache

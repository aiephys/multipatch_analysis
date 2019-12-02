import os
from ..pipeline_module import DatabasePipelineModule


class MultipatchPipelineModule(DatabasePipelineModule):
    
    def job_status(self):
        """Extends DatabasePipelineModule to provide more information about the source of each job.
        """
        db = self.database
        jobs = DatabasePipelineModule.job_status(self)
        expts = {rec[0]:rec[1] for rec in db.query(db.Experiment.ext_id, db.Experiment).all()}
        
        for jid, (status, error, meta) in jobs.items():
            if status is True:
                continue
            if meta is None:
                if jid not in expts:
                    continue
                meta = {'source': expts[jid].original_path}
            else:
                # a path should already be present in meta; translate backward to original source
                sync_file = os.path.join(meta['source'], 'sync_source')
                if os.path.exists(sync_file):
                    meta['source'] = open(sync_file).read()
                    
            jobs[jid] = (status, error, meta)
        
        return jobs

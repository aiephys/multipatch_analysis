import sys, time, multiprocessing, traceback
from datetime import datetime
import numpy as np
from collections import OrderedDict
from pyqtgraph import toposort
from .. import database as db


class PipelineModule(object):
    """Pipeline modules represent analysis tasks that can be run independently of other parts of the analysis
    pipeline. 
    
    For any given experiment, a sequence of analysis stages must be processed in order. Each stage requires a specific
    set of inputs to be present, and produces/stores some output. Inputs and outputs can be raw data files, database tables,
    etc., although outputs are probably always written into a database. Each PipelineModule subclass represents a single stage
    in the sequence of analyses that occur across a pipeline.
    
    The work done by a single stage in the pipeline is divided up into jobs (units of work), where each stage may
    decide for itself what a suitable unit of work is. For many stages, the unit of work is the analysis done on a single experiment.
    Some stages may have finer granularity (multiple jobs per experiment), and some stages may have coarser granularity
    (multiple experiments per job), or even only have a single unit of work for the entire database, such as when aggregating 
    very high-level results.
    
    Note that PipelineModule classes generally need not implement any actual _analysis_; rather, they are simply responsible for
    data _management_ within the analysis pipeline. When an PipelineModule is asked to update an analysis result, it may call out to
    other packages to do the real work.
    
    From here, we should be able to:
    - Ask for which jobs this analysis has already been run (and when)
    - Ask for which jobs this analysis needs to be run (input dependencies are met and no result exists yet or is out of date)    
    - Run and store analysis for any specific experiment
    - Remove / re-run analysis for any specific experiment
    - Remove / re-run analysis for all jobs (for example, after a code change affecting all analysis results)
    - Initialize storage for analysis results (create DB tables, empty folders, etc.)
    
    """    
    
    name = None
    dependencies = []

    @staticmethod
    def all_modules():
        """Return an ordered dictionary mapping {name: module} of all known pipeline modules,
        topologically sorted by dependencies (least dependent to most dependent).
        """
        subclasses = [PipelineModule]
        for sc in subclasses:
            subclasses.extend(sc.__subclasses__())
        
        excluded = [PipelineModule, DatabasePipelineModule]
        deps = {c:c.dependencies for c in subclasses if c not in excluded}
        return OrderedDict([(mod.name, mod) for mod in toposort(deps)])
    
    @classmethod
    def dependent_modules(cls):
        """Return a list of other modules that directly depend on this module.
        """
        mods = cls.all_modules()
        return [mod for mod in mods.values() if cls in mod.dependencies]
    
    @classmethod
    def all_dependent_modules(cls):
        """Return a list of other modules that directly or indirectly depend on this module.
        """
        deps = cls.dependent_modules()
        for dep in deps:
            deps.extend(dep.dependent_modules())
        return [mod for mod in PipelineModule.all_modules().values() if mod in deps]
    
    @classmethod
    def update(cls, job_ids=None, limit=None, parallel=False, workers=None, raise_exceptions=False):
        """Update analysis results for this module.
        
        Parameters
        ----------
        job_ids : list | None
            List of job IDs to be updated, or None to update all jobs.
        parallel : bool
            If True, run jobs in parallel threads or subprocesses.
        workers : int or None
            Number of parallel workers to spawn. If None, then use one worker per CPU core.
        limit : int | None
            Maximum number of jobs to process (or None to disable this limit).
            If limit is enabled, then jobs are randomly shuffled before selecting the limited subset.
        raise_exceptions : bool
            If True, then exceptions are raised and will end any further processing.
            If False, then errors are logged and ignored.
            This is used mainly for debugging to allow traceback inspection.
        """
        print("Updating pipeline stage: %s" % cls.name)
        if job_ids is None:
            print("Searching for jobs to update..")
            run_job_ids, drop_job_ids = cls.updatable_jobs()                    
            if limit is not None:
                np.random.shuffle(job_ids)
                job_ids = job_ids[:limit]
        else:
            run_job_ids = job_ids
            drop_job_ids = job_ids
            
        print("Found %d jobs to update." % len(run_job_ids))

        if len(drop_job_ids) > 0:
            print("Dropping %d invalid results.." % len(drop_job_ids))
            cls.drop_jobs(drop_job_ids)

        run_jobs = [(expt_id, i, len(run_job_ids)) for i, expt_id in enumerate(run_job_ids)]

        if parallel:
            print("Processing all jobs (parallel)..")
            pool = multiprocessing.Pool(processes=workers)
            parallel_jobs = [(cls, job) for job in run_jobs]
            pool.map(run_job_parallel, parallel_jobs)
        else:
            print("Processing all jobs (serial)..")
            for job in run_jobs:
                cls._run_job(job, raise_exceptions=raise_exceptions)

    @classmethod
    def _run_job(cls, job, raise_exceptions=False):
        """Entry point for running a single analysis job; may be invoked in a subprocess.
        """
        job_id, job_index, n_jobs = job
        print("Processing %d/%d  %0.3f") % (job_index, n_jobs, job_id)
        start = time.time()
        try:
            cls.process_job(job_id)
        except Exception:
            if raise_exceptions:
                raise
            else:
                print("Error processing %d/%d  %0.3f:") % (job_index, n_jobs, job_id)
                sys.excepthook(*sys.exc_info())
        else:
            print("Finished %d/%d  %0.3f  (%0.2g sec)") % (job_index, n_jobs, job_id, time.time()-start)
    
    @classmethod
    def process_job(cls, job_id):
        """Process analysis for one job.
        
        Parameters
        ----------
        job_id :
            The ID of the job to be processed.
        
        Must be reimplemented in subclasses.
        """
        raise NotImplementedError()

    @classmethod
    def initialize(cls):
        """Create space (folders, tables, etc.) for this analyzer to store its results.
        """
        raise NotImplementedError()
        
    @classmethod
    def drop_jobs(cls, job_ids):
        """Remove all results previously stored for a list of job IDs.
        """
        raise NotImplementedError()

    @classmethod
    def drop_all(cls):
        """Remove all results generated by this module.
        """
        raise NotImplementedError()

    @classmethod
    def finished_jobs(cls):
        """Return an ordered dict of job IDs that have been processed by this module and
        the dates when they were processed.

        Note that some results returned may be obsolete if dependencies have changed.
        """
        raise NotImplementedError()

    @classmethod
    def ready_jobs(cls):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        jobs = OrderedDict()
        for i,dep in enumerate(cls.dependencies):
            dep_finished = dep.finished_jobs()
            for k,v in dep_finished.items():
                jobs.setdefault(k, []).append(v)
        
        finished = OrderedDict()
        for job, dates in jobs.items():
            if len(dates) < len(cls.dependencies):
                # not all dependencies are finished yet
                continue
            finished[job] = max(dates)

        return finished

    @classmethod
    def updatable_jobs(cls):
        """Return lists of jobs that should be updated and/or should have their results dropped.
        """
        run_job_ids = []
        drop_job_ids = []
        ready = cls.ready_jobs()
        finished = cls.finished_jobs()
        for job in ready:
            if job in finished:
                if ready[job] > finished[job]:
                    # result is invalid
                    run_job_ids.append(job)
                    drop_job_ids
                else:
                    # result is valid
                    pass  
            else:
                # no current result
                run_job_ids.append(job)
        
        # look for orphaned results
        for job in finished:
            if job not in ready and job not in drop_job_ids:
                drop_job_ids.append(job)
        
        print("%d jobs ready for processing, %d finished, %d need drop, %d need update" % (len(ready), len(finished), len(drop_job_ids), len(run_job_ids)))
        return run_job_ids, drop_job_ids


def run_job_parallel(job):
    # multiprocessing Pool.map doesn't work on methods; must be a plain function
    cls, job = job
    cls._run_job(job)


class DatabasePipelineModule(PipelineModule):
    """PipelineModule that implements default behaviors for interacting with database.
    
    * Manages adding / removing entries from pipeline table to indicate job status
    * Default implementations for dropping records 
    """
    table_group = None

    @classmethod
    def create_db_entries(cls, job_id, session):
        """Generate DB entries for *job_id* and add them to *session*.
        """
        raise NotImplementedError()
        
    @classmethod
    def job_query(cls, job_ids, session):
        """Return a query that returns records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        raise NotImplementedError()

    @classmethod
    def dependent_job_ids(cls, module, job_ids):
        """Return a list of all finished job IDs in this module that depend on 
        specific jobs from another module.
        """
        raise NotImplementedError()

    @classmethod
    def process_job(cls, job_id):
        session = db.Session(readonly=False)
        # drop old pipeline job record
        session.query(db.Pipeline).filter(db.Pipeline.job_id==job_id).filter(db.Pipeline.module_name==cls.name).delete()
        session.commit()
        
        try:
            errors = cls.create_db_entries(job_id, session)
            job_result = db.Pipeline(module_name=cls.name, job_id=job_id, success=True, error=errors, finish_time=datetime.now())
            session.add(job_result)

            session.commit()
        except Exception:
            session.rollback()
            
            err = traceback.format_exception(*sys.exc_info())
            job_result = db.Pipeline(module_name=cls.name, job_id=job_id, success=False, error=err, finish_time=datetime.now())
            session.add(job_result)
            session.commit()
            raise
        finally:
            session.close()

    @classmethod
    def initialize(cls):
        """Create space (folders, tables, etc.) for this analyzer to store its results.
        """
        cls.table_group.create_tables()
        
    @classmethod
    def drop_all(cls, reinitialize=True):
        """Remove all results generated by this module.
        """
        for dep in reversed(cls.all_dependent_modules()):
            dep.drop_all(reinitialize=False)
            
        # drop tables and pipeline job records
        cls.table_group.drop_tables()
        session = db.Session(readonly=False)
        session.query(db.Pipeline).filter(db.Pipeline.module_name==cls.name).delete()
        session.commit()

        if reinitialize:
            cls.initialize()        
            for dep in cls.dependent_modules():
                dep.initialize()
    
    @classmethod
    def drop_jobs(cls, job_ids, session=None):
        """Remove all results previously stored for a list of job IDs.
        
        The associated results of dependent modules are also removed.
        """
        if session is None:
            session = db.Session(readonly=False)
        
        for dep in reversed(cls.dependent_modules()):
            dep_jobs = dep.dependent_job_ids(cls, job_ids)
            dep.drop_jobs(dep_jobs, session=session)
        
        jobs = cls.job_query(job_ids, session).all()
        for job in jobs:
            session.delete(job)
        session.query(db.Pipeline).filter(db.Pipeline.module_name==cls.name).filter(db.Pipeline.job_id.in_(job_ids)).delete(synchronize_session=False)
        session.commit()
        
        print("Dropped %d jobs from %s module" % (len(jobs), cls.name))
    
    @classmethod
    def finished_jobs(cls):
        """Return an ordered dict of job IDs that have been processed by this module and
        the dates when they were processed.

        Note that some results returned may be obsolete if dependencies have changed.
        """
        session = db.Session()
        slices = session.query(db.Pipeline.job_id, db.Pipeline.finish_time).filter(db.Pipeline.module_name==cls.name).filter(db.Pipeline.success==True).all()
        session.rollback()
        return OrderedDict([(uid, date) for uid, date in slices])

import os, glob, re, pickle, time
from datetime import datetime
from collections import OrderedDict
from acq4.util.DataManager import getDirHandle
from ..pipeline_module import DatabasePipelineModule
from ... import config, lims, constants
from ...util import datetime_to_timestamp, timestamp_to_datetime
from ...data.slice import Slice


class SlicePipelineModule(DatabasePipelineModule):
    """Imports per-slice metadata into DB.
    """
    name = 'slice'
    dependencies = []
    table_group = ['slice']
    
    @classmethod
    def create_db_entries(cls, job, session):
        job_id = job['job_id']
        db = job['database']

        slices = all_slices()
        path = slices[job_id]
        
        sl = Slice.get(path)

        fields = {
            'ext_id': job_id,
            'acq_timestamp': sl.timestamp,
            'species': sl.species,
            'date_of_birth': sl.date_of_birth,
            'age': sl.age,
            'sex': sl.sex,
            'genotype': None if sl.genotype is None else sl.genotype.gtype,
            'orientation': sl.orientation,
            'surface': sl.surface,
            'hemisphere': sl.hemisphere,
            'quality': sl.quality,
            'slice_time': sl.slice_time,
            'slice_conditions': {},
            'lims_specimen_name': sl.lims_specimen_name,
            'storage_path': sl.storage_path,
        }

        sl = db.Slice(**fields)
        session.add(sl)

    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        """
        db = self.database
        return session.query(db.Slice).filter(db.Slice.acq_timestamp.in_(job_ids)).all()

    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        slices = all_slices()
        ready = OrderedDict()
        for ts, path in slices.items():
            mtime = os.stat(os.path.join(path, '.index')).st_mtime
            # test file updates:
            # import random
            # if random.random() > 0.8:
            #     mtime *= 2
            ready[ts] = timestamp_to_datetime(mtime)
        return ready


_all_slices = None
def all_slices():
    """Return a dict mapping {slice_timestamp: path} for all known slices.
    
    This is only generated once per running process; set _all_slices = None
    to force the list to be regenerated.
    """
    global _all_slices
    if _all_slices is not None:
        return _all_slices
        
    # Speed things up by caching this list with a 4 hour timeout
    cachefile = os.path.join(config.cache_path, 'all_slices.pkl')
    if os.path.exists(cachefile):
        age = time.time() - os.stat(cachefile).st_mtime
        if age < 4 * 3600:
            print("Loaded slice timestamps from cache (%0.1f hours old)" % (age/3600.))
            return pickle.load(open(cachefile, 'rb'))
    
    slice_dirs = sorted(glob.glob(os.path.join(config.synphys_data, '*', 'slice_*')))

    _all_slices = OrderedDict()
    for path in slice_dirs:
        dh = getDirHandle(path)
        ts = dh.info().get('__timestamp__')
        if ts is None:
            print("MISSING TIMESTAMP: %s" % path)
            continue
        _all_slices["%0.3f"%ts] = path
        
    try:
        tmpfile = cachefile+'.tmp'
        pickle.dump(_all_slices, open(tmpfile, 'wb'))
        os.rename(tmpfile, cachefile)
    except:
        if os.path.exists(tmpfile):
            os.remove(tmpfile)
    
    return _all_slices

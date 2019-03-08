import os, glob, re
from datetime import datetime
from collections import OrderedDict
from acq4.util.DataManager import getDirHandle
from .. import database as db
from ..database import slice_tables
from .pipeline_module import DatabasePipelineModule
from .. import config
from .. import lims
from ..util import datetime_to_timestamp


class SlicePipelineModule(DatabasePipelineModule):
    """Imports per-slice metadata into DB.
    """
    name = 'slice'
    dependencies = []
    table_group = slice_tables
    
    @classmethod
    def process_job(cls, job_id):
        slices = all_slices()
        path = slices[job_id]
        dh = getDirHandle(path)
        info = dh.info()
        parent_info = dh.parent().info()
        
        # pull some metadata from LIMS
        sid = info['specimen_ID'].strip()
        limsdata = lims.specimen_info(sid)

        quality = info.get('slice quality', None)
        try:
            quality = int(quality)
        except Exception:
            quality = None

        # Interpret slice time
        slice_time = parent_info.get('time_of_dissection', None)
        if slice_time is not None:
            m = re.match(r'((20\d\d)-(\d{1,2})-(\d{1,2}) )?(\d+):(\d+)', slice_time.strip())
            if m is not None:
                _, year, mon, day, hh, mm = m.groups()
                if year is None:
                    date = datetime.fromtimestamp(dh.parent().info()['__timestamp__'])
                    slice_time = datetime(date.year, date.month, date.day, int(hh), int(mm))
                else:
                    slice_time = datetime(int(year), int(mon), int(day), int(hh), int(mm))

        fields = {
            'acq_timestamp': info['__timestamp__'],
            'species': limsdata['organism'],
            'date_of_birth': limsdata['date_of_birth'],
            'age': limsdata['age'],
            'sex': limsdata['sex'],
            'genotype': limsdata['genotype'],
            'orientation': limsdata['plane_of_section'],
            'surface': limsdata['exposed_surface'],
            'hemisphere': limsdata['hemisphere'],
            'quality': quality,
            'slice_time': slice_time,
            'slice_conditions': {},
            'lims_specimen_name': sid,
            'storage_path': dh.name(relativeTo=dh.parent().parent()),
        }

        sl = db.Slice(**fields)
        # sl.meta = {'db_timestamp': time.time()}
        
        session = db.Session(readonly=False)
        try:
            session.add(sl)
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
        slices = session.query(db.Slice).filter(db.Slice.acq_timestamp.in_(job_ids))
        for sl in slices:
            session.delete(sl)
        session.commit()

    @classmethod
    def finished_jobs(cls):
        """Return an ordered dict of job IDs that have been processed by this module and
        the dates when they were processed.

        Note that some results returned may be obsolete if dependencies have changed.
        """
        session = db.Session()
        slices = session.query(db.Slice.acq_timestamp, db.Slice.time_created).all()
        session.rollback()
        return OrderedDict([(uid, datetime_to_timestamp(date)) for uid, date in slices])

    @classmethod
    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        slices = all_slices()
        ready = OrderedDict()
        for ts, path in slices.items():
            age = os.stat(os.path.join(path, '.index')).st_mtime
            ready[ts] = age
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
    
    slice_dirs = sorted(glob.glob(os.path.join(config.synphys_data, '15034*', 'slice_*')))
    _all_slices = OrderedDict()
    for path in slice_dirs:
        dh = getDirHandle(path)
        ts = dh.info()['__timestamp__']
        _all_slices[ts] = path
        
    return _all_slices

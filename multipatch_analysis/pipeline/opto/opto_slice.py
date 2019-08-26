from multipatch_analysis.pipeline.pipeline_module import DatabasePipelineModule
from optoanalysis import data_model
#from multipatch_analysis.database import slice_tables
import os, glob, re, pickle, time, csv
from multipatch_analysis import config, lims, constants
from acq4.util.DataManager import getDirHandle
from collections import OrderedDict
from multipatch_analysis.util import datetime_to_timestamp, timestamp_to_datetime
#import multipatch_analysis.database as db


#from .opto_experiment import OptoExperimentPipelineModule

class OptoSlicePipelineModule(DatabasePipelineModule):

    name = 'opto_slice'
    depencencies = []
    #table_group = slice_tables
    table_group = ['slice']

    @classmethod
    def create_db_entries(cls, job, session):
        job_id = job['job_id']
        db = job['database']

        slices = all_slices()
        path = slices[job_id]

        if path == 'place_holder':
            sl = db.Slice(storage_path='place_holder', acq_timestamp=0.0)
            session.add(sl)
            return

        dh = getDirHandle(path)
        info = dh.info()
        parent_info = dh.parent().info()
        
        # pull some metadata from LIMS
        #sid = self.find_specimen_name(dh)
        sids = data_model.find_lims_specimen_ids(dh)
        #print('sids:', sids)
        if len(sids) == 0:
            limsdata = {}
        elif len(sids) == 1:
            limsdata = lims.specimen_info(specimen_id=sids[0])
        elif len(sids) > 1:
            data = []
            for i in sids:
                data.append(lims.specimen_info(specimen_id=i))
            limsdata = {}
            for key in ['organism', 'date_of_birth', 'age', 'sex', 'plane_of_section', 'exposed_surface', 'hemisphere', 'specimen_name', 'genotype']:
                vals = set([d[key] for d in data])
                if len(vals) == 1:
                    limsdata[key] = vals[0]


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

        # construct full genotype string 
        genotype = limsdata.get('genotype', '')
        for info in (parent_info, info):
            inj = info.get('injections')
            if inj in (None, ''):
                continue
            if inj not in constants.INJECTIONS:
                raise KeyError("Injection %r is unknown in constants.INJECTIONS" % inj)
            genotype = genotype + ';' + constants.INJECTIONS[inj]


        fields = {
            'acq_timestamp': info['__timestamp__'],
            'species': limsdata.get('organism'),
            'date_of_birth': limsdata.get('date_of_birth'),
            'age': limsdata.get('age'),
            'sex': limsdata.get('sex'),
            'genotype': genotype,
            'orientation': limsdata.get('plane_of_section'),
            'surface': limsdata.get('exposed_surface'),
            'hemisphere': limsdata.get('hemisphere'),
            'quality': quality,
            'slice_time': slice_time,
            'slice_conditions': {},
            'lims_specimen_name': limsdata.get('specimen_name'),
            'storage_path': dh.name(relativeTo=getDirHandle(config.synphys_data)),
        }

        sl = db.Slice(**fields)
        session.add(sl)
        session.commit()

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
            if path == 'place_holder':
                ready['%.3f'%0.0] = os.stat(config.experiment_csv).st_mtime
            else:
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
            return pickle.load(open(cachefile, 'r'))
    
    #slice_dirs = sorted(glob.glob(os.path.join(config.synphys_data, '*', 'slice_*')))
    expt_csv = config.experiment_csv
    csv_entries = []
    with open(expt_csv, 'r') as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            csv_entries.append(row)

    slice_dirs = sorted([os.path.split(os.path.join(config.synphys_data, exp['site_path']))[0] for exp in csv_entries])

    _all_slices = OrderedDict()
    for path in slice_dirs:
        dh = getDirHandle(path)
        ts = dh.info().get('__timestamp__')
        if ts is None:
            #print("MISSING TIMESTAMP: %s" % path)
            _all_slices.update([('%.3f'%0.0, 'place_holder')])
            continue
        ts = '%0.3f'%ts ## convert timestamp to string here, make sure it has 3 decimal places
        _all_slices[ts] = path
        
    try:
        tmpfile = cachefile+'.tmp'
        pickle.dump(_all_slices, open(tmpfile, 'w'))
        os.rename(tmpfile, cachefile)
    except:
        if os.path.exists(tmpfile):
            os.remove(tmpfile)
    
    return _all_slices

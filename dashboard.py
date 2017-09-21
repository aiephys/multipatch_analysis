import os, datetime, re
import config
import database
from acq4.util.DataManager import getDirHandle



for path in config.raw_data_paths:
    for fname in os.listdir(path):
        fname = os.path.join(path, fname)
        if not os.path.isdir(fname):
            continue
        
        dh = getDirHandle(fname)
        if not dh.isManaged():
            continue
        
        for sl in dh.ls():
            if not sl.startswith('slice_'):
                continue
            slice_dh = dh[sl]
            if not slice_dh.isDir():
                continue
            
            for site in slice_dh.ls():
                if not site.startswith('site_'):
                    continue
                site_dh = slice_dh[site]
                if not site_dh.isDir():
                    continue
                
                ts = site_dh.info()['__timestamp__']
                date = datetime.datetime.fromtimestamp(ts)
                site_id = date
                
                rig = re.search('(MP\d)_', site_dh.name()).groups()[0]
                
                try:
                    expt_entry = database.experiment_from_timestamp(date)
                    expt_steps = expt_entry.submission_data
                    if expt_steps is None:
                        expt_steps = {}
                except KeyError:
                    expt_entry = None
                    expt_steps = {}
                
                print("{rig} {date} {uid} {in_db} {in_server}".format(
                    rig=rig,
                    date=date.strftime('%Y-%m-%d'), 
                    uid='%0.2f'%ts,
                    in_db=expt_entry is not None,
                    in_server='raw_data_location' in expt_steps,
                ))

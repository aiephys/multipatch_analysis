"""
Used to synchronize raw data from rigs to central server.
"""

import os, shutil, glob, traceback, pickle
from acq4.util.DataManager import getDirHandle

import config


class RawDataSubmission(object):
    """Copies all raw data to a central server.
    """
    message = "Copying data to server"
    
    def __init__(self, site_dh):
        self.changes = None
        self.site_dh = site_dh
        
    def check(self):
        return [], []
        
    def summary(self):
        return {'raw_data': None}
        
    def submit(self):
        self.changes = []
        site_dh = self.site_dh
        slice_dh = site_dh.parent()
        expt_dh = slice_dh.parent()
        
        print("Copying %s to server.." % site_dh.name())
        
        # Decide how the top-level directory will be named on the remote server
        # (it may already be there from a previous slice/site, or the current
        # name may already be taken by another rig.)
        server_expt_path = get_experiment_server_path(expt_dh)
        
        self.server_path = server_expt_path
        print("   using server path: %s" % server_expt_path)
        self._sync_paths(expt_dh.name(), server_expt_path)
        
        # Copy slice files if needed
        server_slice_path = os.path.join(server_expt_path, slice_dh.shortName())
        self._sync_paths(slice_dh.name(), server_slice_path)

        # Copy site files if needed
        server_site_path = os.path.join(server_slice_path, site_dh.shortName())
        self._sync_paths(site_dh.name(), server_site_path)
        
        print("Done.")
        
    def _sync_paths(self, source, target):
        """Non-recursive directory sync.
        """
        if not os.path.isdir(target):
            os.mkdir(target)
            self.changes.append(('mkdir', source, target))
        for fname in os.listdir(source):
            src_path = os.path.join(source, fname)
            if os.path.isfile(src_path):
                dst_path = os.path.join(target, fname)
                # don't copy again if destination file is newer than source file
                if os.path.isfile(dst_path):
                    src_mtime = os.stat(src_path).st_mtime
                    dst_mtime = os.stat(dst_path).st_mtime
                    if dst_mtime >= src_mtime:
                        print("    skip %s => %s" % (src_path, dst_path))
                        continue
                    print("    updt %s => %s" % (src_path, dst_path))
                    safe_copy(src_path, dst_path)
                    self.changes.append(('update', src_path, dst_path))
                else:
                    print("    copy %s => %s" % (src_path, dst_path))
                    safe_copy(src_path, dst_path)
                    self.changes.append(('copy', src_path, dst_path))


def get_experiment_server_path(dh):
    server_path = config.synphys_data
    acq_timestamp = dh.info()['__timestamp__']
    
    # First check the cache
    cache = experiment_path_cache()
    if acq_timestamp in cache:
        return cache[acq_timestamp]
    
    # We have not already submitted a site from this experiment folder;
    # look for a suitable new directory name on the server
    expt_base_name = dh.shortName().split('_')[0]
    expt_dirs = set(os.listdir(server_path))
    i = 0
    while True:
        expt_name = expt_base_name + '_%03d' % i
        if expt_name not in expt_dirs:
            break
        i += 1
    server_expt_path = os.path.join(server_path, expt_name)
    assert not os.path.exists(server_expt_path)
    
    os.mkdir(server_expt_path)
    try:
        dh = getDirHandle(server_expt_path)
        dh.setInfo(__timestamp__=acq_timestamp)
    except Exception:
        if os.path.exists(server_expt_path):
            shutil.rmtree(server_expt_path)
        raise
    
    cache[acq_timestamp] = server_expt_path
    write_expt_path_cache()
    
    return server_expt_path


_expt_path_cache = None
def experiment_path_cache():
    global _expt_path_cache
    if _expt_path_cache is None:
        cache_file = os.path.join(config.synphys_data, 'experiment_path_cache.pkl')
        if os.path.isfile(cache_file):
            try:
                _expt_path_cache = pickle.load(open(cache_file, 'rb'))
            except Exception:
                print("Error loading experiment path cache; will regenerate:")
                sys.excepthook(*sys.exc_info())
                generate_expt_path_cache()
        else:
            generate_expt_path_cache()
    return _expt_path_cache


def generate_expt_path_cache():
    global _expt_path_cache
    _expt_path_cache = {}
    root = getDirHandle(config.synphys_data)
    for f in root.ls():
        dh = root[f]
        if not dh.isDir():
            continue
        acq_timestamp = dh.info()['__timestamp__']
        if acq_timestamp in _expt_path_cache:
            raise Exception("timestamp %s appears twice in synphys data!!" % acq_timestamp)
        _expt_path_cache[acq_timestamp] = dh.name()
    write_expt_path_cache()
    
    
def write_expt_path_cache():
    global _expt_path_cache
    cache_file = os.path.join(config.synphys_data, 'experiment_path_cache.pkl')
    tmp = cache_file+'.tmp'
    pickle.dump(_expt_path_cache, open(tmp, 'wb'))
    os.rename(tmp, cache_file)


def safe_copy(src, dst):
    tmp_dst = dst + '_uploading'
    try:
        shutil.copyfile(src, tmp_dst)
        os.rename(tmp_dst, dst)
    finally:
        if os.path.isfile(tmp_dst):
            os.remove(tmp_dst)
    

def sync_experiment(site_dir):
    dh = getDirHandle(site_dir)
    sub = RawDataSubmission(dh)
    err, warn = sub.check()
    if len(err) > 0:
        return [], err, warn
    sub.submit()
    return sub.changes, err, warn
    
    
def find_all_sites(root):
    sites = glob.glob(os.path.join(root, '*', 'slice_*', 'site_*'))
    sites.sort(reverse=True)
    return sites
    

def sync_all():
    log = []
    for raw_data_path in config.raw_data_paths:
        for site_dir in find_all_sites(raw_data_path):
            try:
                changes, err, warn = sync_experiment(site_dir)
                if len(changes) > 0:
                    log.append((site_dir, changes, err, warn))
            except Exception:
                exc = traceback.format_exception()
                log.append((site_dir, [], exc, []))
                



if __name__ == '__main__':
    sync_all()
    
    #sites = find_all_sites(config.raw_data_paths[1])
    #path = get_experiment_server_path(getDirHandle(sites[0]).parent().parent())
    #changes, err, warn = sync_experiment(sites[0])
    
    #print("CHANGES:")
    #for ch in changes:
        #print(ch)
        
    #print("ERRORS:")
    #print("\n".join(err))
    
    #print("WARNINGS:")
    #print("\n".join(warn))
    
    
"""
Used to synchronize raw data from rigs to central server.
"""

import os, sys, shutil, glob, traceback, pickle, time
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
        
        now = time.strftime('%Y-%m-%d_%H:%M:%S')
        self.log("========== %s : Sync %s to server" % (now, site_dh.name()))
        self.skipped = 0
        
        try:
            # Decide how the top-level directory will be named on the remote server
            # (it may already be there from a previous slice/site, or the current
            # name may already be taken by another rig.)
            server_expt_path = get_experiment_server_path(expt_dh)
            
            self.server_path = server_expt_path
            self.log("    using server path: %s" % server_expt_path)
            self._sync_paths(expt_dh.name(), server_expt_path)
            
            # Copy slice files if needed
            server_slice_path = os.path.join(server_expt_path, slice_dh.shortName())
            self._sync_paths(slice_dh.name(), server_slice_path)

            # Copy site files if needed
            server_site_path = os.path.join(server_slice_path, site_dh.shortName())
            self._sync_paths(site_dh.name(), server_site_path)
            
            self.log("    Done; skipped %d files." % self.skipped)
            
            # Leave a note about the source of this data
            open(os.path.join(server_site_path, 'sync_source'), 'wb').write(site_dh.name())
        except Exception:
            err = traceback.format_exc()
            self.changes.append(('error', site_dh.name(), err))
            self.log(err)

    def log(self, msg):
        print(msg)
        with open(os.path.join(config.synphys_data, 'sync_log'), 'ab') as log_fh:
            log_fh.write(msg+'\n')
        
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
                
                # Skip large files:
                #   - pxp > 10GB
                #   - others > 5GB
                src_stat = os.stat(src_path)
                if (src_stat.st_size > 5e9 and not src_path.endswith('.pxp')) or  (src_stat.st_size > 15e9):
                    self.log("    err! %s => %s" % (src_path, dst_path))
                    self.changes.append(('error', src_path, 'file too large'))
                    continue
                
                status = sync_file(src_path, dst_path)
                if status == 'skip':
                    self.skipped += 1
                elif status == 'copy':
                    self.log("    copy %s => %s" % (src_path, dst_path))
                    self.changes.append(('copy', src_path, dst_path))
                elif status == 'update':
                    self.log("    updt %s => %s" % (src_path, dst_path))
                    self.changes.append(('update', src_path, dst_path))


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
        # temporarily mark with timestamp; should be overwritten later.
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


def sync_file(src, dst):
    """Safely copy *src* to *dst*, but only if *src* is newer or a different size.
    """
    if os.path.isfile(dst):
        src_stat = os.stat(src)
        dst_stat = os.stat(dst)
        up_to_date = dst_stat.st_mtime >= src_stat.st_mtime and src_stat.st_size == dst_stat.st_size
        
        if up_to_date:
            return "skip"
        
        safe_copy(src, dst)
        return "update"
    else:
        safe_copy(src, dst)
        return "copy"


def safe_copy(src, dst):
    """Copy a file, but rename the destination file if it already exists.
    
    Also, the destination file is suffixed ".partial" until the copy is complete.
    """
    tmp_dst = dst + '.partial'
    try:
        new_name = None
        chunk_copy(src, tmp_dst)
        if os.path.exists(dst):
            # rename destination file to avoid overwriting
            now = time.strftime('%Y-%m-%d_%H:%M:%S')
            i = 0
            while True:
                new_name = '%s_%s_%d' % (dst, now, i)
                if not os.path.exists(new_name):
                    break
                i += 1
            os.rename(dst, new_name)
        os.rename(tmp_dst, dst)
    except Exception:
        # Move dst file back if there was a problem during copy
        if new_name is not None and os.path.exists(new_name):
            os.rename(new_name, dst)
        raise
    finally:
        if os.path.isfile(tmp_dst):
            os.remove(tmp_dst)
    

def chunk_copy(src, dst, chunk_size=100e6):
    """Manually copy a file one chunk at a time.
    
    This allows progress feedback and more graceful cancellation during long
    copy operations.
    """
    if os.path.exists(dst):
        raise Exception("Won't copy over existing file %s" % dst)
    size = os.stat(src).st_size
    in_fh = open(src, 'rb')
    out_fh = open(dst, 'ab')
    msglen = 0
    try:
        with in_fh:
            with out_fh:
                chunk_size = int(chunk_size)
                tot = 0
                while True:
                    chunk = in_fh.read(chunk_size)
                    out_fh.write(chunk)
                    tot += len(chunk)
                    if size > chunk_size * 2:
                        n = int(50 * (float(tot) / size))
                        msg = ('[' + '#' * n + '-' * (50-n) + ']  %d / %d MB\r') % (int(tot/1e6), int(size/1e6))
                        msglen = len(msg)
                        sys.stdout.write(msg)
                        sys.stdout.flush()
                    if len(chunk) < chunk_size:
                        break
                sys.stdout.write("[###  flushing..  \r")
                sys.stdout.flush()
        sys.stdout.write(' '*msglen + '\r')
        sys.stdout.flush()
    except Exception:
        if os.path.isfile(dst):
            os.remove(dst)
        raise
    
    

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
    

def sync_all(log):
    for raw_data_path in config.raw_data_paths:
        paths = find_all_sites(raw_data_path)
        sync_paths(paths, log)


def sync_paths(paths, log):
    for site_dir in paths:
        try:
            changes, err, warn = sync_experiment(site_dir)
            if len(changes) > 0:
                log.append((site_dir, changes, err, warn))
        except Exception:
            exc = traceback.format_exc()
            print(exc)
            log.append((site_dir, [], exc, []))


if __name__ == '__main__':
    log = []
    
    paths = sys.argv[1:]
    if len(paths) == 0:
        sync_all(log)
    else:
        sync_paths(paths, log)
    
    errs = [change for site in log for change in site[1] if change[0] == 'error']
    print("\n----- DONE ------\n   %d errors" % len(errs))
    
    for err in errs:
        print(err[1], err[2])
    
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
    
    
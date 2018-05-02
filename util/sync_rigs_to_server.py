"""
Used to synchronize raw data from rigs to central server.


Practical note: The file synchronization implements an rsync-like
functionality in Python. In theory we would have preferred to just use 
rsync itself, but in practice this presents a few difficult issues:

* Rsync on windows requires cygwin and sshd. These work in most cases,
  but they are not as thoroughly developed as their linux counterparts.
  Bugs encountered in this process can be showstoppers.
* The actual file synchronization code here is simple and transparent,
  allowing us to customize versioning, error checking, logging,
  and file name transformations. Achieving the equivalent with rsync
  actually requires more code than we have written here, just to handle
  subprocessing, CLI flag generation, and fragile pipe communication.
"""

import os, sys, shutil, glob, traceback, pickle, time
from acq4.util.DataManager import getDirHandle

from multipatch_analysis import config
from multipatch_analysis.util import sync_file


def sync_experiment(site_dir):
    """Synchronize all files for an experiment to the server.

    Return a list of changes made.
    """
    site_dh = getDirHandle(site_dir)
    changes = []
    slice_dh = site_dh.parent()
    expt_dh = slice_dh.parent()
    
    now = time.strftime('%Y-%m-%d_%H:%M:%S')
    log("========== %s : Sync %s to server" % (now, site_dh.name()))
    skipped = 0
    
    try:
        # Decide how the top-level directory will be named on the remote server
        # (it may already be there from a previous slice/site, or the current
        # name may already be taken by another rig.)
        server_expt_path = get_experiment_server_path(expt_dh)
        
        log("    using server path: %s" % server_expt_path)
        skipped += _sync_paths(expt_dh.name(), server_expt_path, changes)
        
        # Copy slice files if needed
        server_slice_path = os.path.join(server_expt_path, slice_dh.shortName())
        skipped += _sync_paths(slice_dh.name(), server_slice_path, changes)

        # Copy site files if needed
        server_site_path = os.path.join(server_slice_path, site_dh.shortName())
        skipped += _sync_paths(site_dh.name(), server_site_path, changes)
        
        log("    Done; skipped %d files." % skipped)
        
        # Leave a note about the source of this data
        open(os.path.join(server_site_path, 'sync_source'), 'wb').write(site_dh.name())
    except Exception:
        err = traceback.format_exc()
        changes.append(('error', site_dh.name(), err))
        log(err)

    return changes


def log(msg):
    print(msg)
    with open(os.path.join(config.synphys_data, 'sync_log'), 'ab') as log_fh:
        log_fh.write(msg+'\n')


def _sync_paths(source, target, changes):
    """Non-recursive directory sync.

    Return the number of skipped files.
    """
    skipped = 0
    if not os.path.isdir(target):
        os.mkdir(target)
        changes.append(('mkdir', source, target))
    for fname in os.listdir(source):
        src_path = os.path.join(source, fname)
        if os.path.isfile(src_path):
            dst_path = os.path.join(target, fname)
            
            # Skip large files:
            #   - pxp > 10GB
            #   - others > 5GB
            src_stat = os.stat(src_path)
            if (src_stat.st_size > 5e9 and not src_path.endswith('.pxp')) or  (src_stat.st_size > 15e9):
                log("    err! %s => %s" % (src_path, dst_path))
                changes.append(('error', src_path, 'file too large'))
                continue
            
            status = sync_file(src_path, dst_path)
            if status == 'skip':
                skipped += 1
            elif status == 'copy':
                log("    copy %s => %s" % (src_path, dst_path))
                changes.append(('copy', src_path, dst_path))
            elif status == 'update':
                log("    updt %s => %s" % (src_path, dst_path))
                changes.append(('update', src_path, dst_path))

    return skipped


def get_experiment_server_path(dh):
    server_path = config.synphys_data
    acq_timestamp = dh.info()['__timestamp__']
    
    # First check the cache
    cache = experiment_path_cache()
    if acq_timestamp in cache:
        return os.path.join(server_path, cache[acq_timestamp])
    
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
    
    cache[acq_timestamp] = expt_name
    write_expt_path_cache()
    
    return server_expt_path


_expt_path_cache = None
def experiment_path_cache():
    """The server contains a pickled dictionary mapping expt_uid:path
    """
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
        if 'recycle' in f.lower():
            continue
        dh = root[f]
        if not dh.isDir():
            continue
        try:
            acq_timestamp = dh.info()['__timestamp__']
        except KeyError:
            print("NO TIMESTAMP:", dh.name())
            sys.exit(-1)
        if acq_timestamp in _expt_path_cache:
            raise Exception("timestamp %s appears twice in synphys data!!" % acq_timestamp)
        _expt_path_cache[acq_timestamp] = dh.name(relativeTo=root)
    write_expt_path_cache()
    
    
def write_expt_path_cache():
    global _expt_path_cache
    cache_file = os.path.join(config.synphys_data, 'experiment_path_cache.pkl')
    tmp = cache_file+'.tmp'
    pickle.dump(_expt_path_cache, open(tmp, 'wb'))
    os.rename(tmp, cache_file)
    

def find_all_sites(root):
    sites = glob.glob(os.path.join(root, '*', 'slice_*', 'site_*'))
    sites.sort(reverse=True)
    return sites
    

def sync_all():
    """Synchronize all known rig data paths to the server
    """
    log = []
    synced_paths = []
    # Loop over all rigs
    for rig_name, data_paths in config.rig_data_paths.items():
        # Each rig may have multiple paths to check
        for data_path in data_paths:
            # each "primary" storage path also has a corresponding "archive" path
            # on the same machine, but we probably only need to synchronize from the
            # primary storage.
            data_path = data_path['primary']

            # Get a list of all experiments stored in this path
            paths = find_all_sites(data_path)

            # synchronize files for each experiment to the server
            log.extend(sync_paths(paths))
            synced_paths.append((rig_name, data_path))
    
    return log, synced_paths


def sync_paths(paths):
    log = []
    for site_dir in paths:
        try:
            changes = sync_experiment(site_dir)
            if len(changes) > 0:
                log.append((site_dir, changes))
        except Exception:
            exc = traceback.format_exc()
            print(exc)
            log.append((site_dir, [], exc, []))
    return log


if __name__ == '__main__':
    
    paths = sys.argv[1:]
    if len(paths) == 0:
        # Synchronize all known rig data paths
        log, synced_paths = sync_all()
        print("==========================\nSynchronized files from:")
        for rig_name, data_path in synced_paths:
            print("%s  :  %s" % (rig_name, data_path))

    else:
        # synchronize just the specified path(s)
        log = sync_paths(paths)
    
    errs = [change for site in log for change in site[1] if change[0] == 'error']
    print("\n----- DONE ------\n   %d errors" % len(errs))
    
    for err in errs:
        print(err[1], err[2])

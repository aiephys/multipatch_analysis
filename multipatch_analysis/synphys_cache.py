import os, sys, glob
import config
from .util import sync_file


_cache = None
def get_cache():
    global _cache
    if _cache is None:
        _cache = SynPhysCache()
    return _cache


class SynPhysCache(object):
    """Maintains a local cache of files from the synphys raw data repository.
    """
    def __init__(self, local_path=config.cache_path, remote_path=config.synphys_data):
        # If a relative path is given, then interpret it as relative to home
        if not os.path.isabs(local_path):
            local_path = os.path.join(os.path.dirname(__file__), '..', local_path)

        self.local_path = os.path.abspath(local_path)
        self.remote_path = os.path.abspath(remote_path)
        
    def list_nwbs(self):
        return glob.glob(os.path.join(self.remote_path, '*', 'slice_*', 'site_*', '*.nwb'))
    
    def list_pip_yamls(self):
        return glob.glob(os.path.join(self.remote_path, '*', 'slice_*', 'site_*', 'pipettes.yml'))

    def get_cache(self, filename):
        filename = os.path.abspath(filename)
        if not filename.startswith(self.remote_path):
            raise Exception("Requested file %s is not inside %s" % (filename, self.remote_path))

        rel_filename = filename[len(self.remote_path):].lstrip(os.sep)
        path, _ = os.path.split(rel_filename)
        
        local_path = os.path.join(self.local_path, path)
        self._mkdir(local_path)
        
        local_filename = os.path.join(self.local_path, rel_filename)
        
        sync_file(filename, local_filename)
        return local_filename
        
    def _mkdir(self, path):
        if not os.path.isdir(path):
            root, _ = os.path.split(path)
            if root != '':
                self._mkdir(root)
            os.mkdir(path)

import os, sys, glob
from collections import OrderedDict
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
        
    def list_experiments(self):
        yamls = self.list_pip_yamls()
        site_dirs = sorted([os.path.dirname(yml) for yml in yamls], reverse=True)
        expts = OrderedDict([(dir_timestamp(site_dir), site_dir) for site_dir in site_dirs])
        return expts

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


def dir_timestamp(path):
    """Get the timestamp from an index file.

    This is just a very lightweight version of the same functionality provided by ACQ4's DirHandle.info()['__timestamp__'].
    We'd prefer not to duplicate this functionality, but acq4 has UI dependencies that make automated scripting more difficult.
    """
    index_file = os.path.join(path, '.index')
    in_dir = False
    search_indent = None
    for line in open(index_file, 'rb').readlines():
        if line.startswith('.:'):
            in_dir = True
            continue
        if line[0] != ' ':
            if in_dir is True:
                return None
        if not in_dir:
            continue
        indent = len(line) - len(line.lstrip(' '))
        if search_indent is None:
            search_indent = indent
        if indent != search_indent:
            continue
        line = line.lstrip()
        key = '__timestamp__:'
        if line.startswith(key):
            return float(line[len(key):])

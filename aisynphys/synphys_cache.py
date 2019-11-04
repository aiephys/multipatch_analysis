import os, sys, glob, pickle, base64, urllib
from collections import OrderedDict
from . import config
from .util import sync_file, dir_timestamp, interactive_download


_db_versions = None
def list_db_versions():
    """Return a dictionary listing database versions that are available for download.
    """
    global _db_versions
    if _db_versions is None:
        # DB urls are stored as a base64-encoded pickle on GitHub.
        # This allows us to change download URLs without requiring users to pull new code.
        # The b64 encoding is just intended to prevent bots scraping our URLs
        b64_urls = urllib.request.urlopen('https://raw.githubusercontent.com/AllenInstitute/aisynphys/download_urls/download_urls').read()
        _db_versions = pickle.loads(base64.b64decode(b64_urls))
    return _db_versions


def get_db_path(db_version):
    """Return the filesystem path of a known database file.
    
    If the file does not exist locally, then it will be downloaded before returning
    the path.
    """
    cache_path = os.path.join(config.cache_path, 'database')
    cache_file = os.path.join(cache_path, db_version)
    
    if not os.path.exists(cache_path):
        os.makedirs(cache_path)
    if not os.path.exists(cache_file):
        versions = list_db_versions()
        if db_version not in versions:
            raise KeyError("Unknown database version; options are: %s" % str(list(versions.keys())))
        url = versions[db_version]['url']
        interactive_download(url, cache_file)
        
    return cache_file

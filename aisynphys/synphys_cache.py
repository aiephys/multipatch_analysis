import os, sys, glob, pickle, base64, urllib, json, re
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


_file_index = None
def get_data_file_index():
    global _file_index
    if _file_index is None:
        query_url = "http://api.brain-map.org/api/v2/data/WellKnownFile/query.json?criteria=[path$il*synphys*]&num_rows=%d"

        # request number of downloadable files
        count_json = urllib.request.urlopen(query_url % 0).read()
        count = json.loads(count_json)
        if not count['success']:
            raise Exception("Error loading file index: %s" % count['msg'])

        # request full index
        index_json = urllib.request.urlopen(query_url % count['total_rows']).read()
        index = json.loads(index_json)
        if not index['success']:
            raise Exception("Error loading file index: %s" % index['msg'])

        # extract {expt_id:url} mapping from index
        _file_index = {}
        for rec in index['msg']:
            m = re.match(r'.*-(\d+\.\d+)\.nwb$', rec['path'])
            if m is None:
                # skip non-nwb files
                continue
            expt_id = m.groups()[0]
            _file_index[expt_id] = rec['download_link']

    return _file_index


def get_nwb_path(expt_id):
    """Return the local filesystem path to an experiment's nwb file. 

    If the file does not exist locally, then attempt to download.
    """
    cache_path = os.path.join(config.cache_path, 'raw_data_files', expt_id)
    cache_file = os.path.join(cache_path, 'data.nwb')
    
    if not os.path.exists(cache_path):
        os.makedirs(cache_path)
    if not os.path.exists(cache_file):
        index = get_data_file_index()
        url = index.get(expt_id, None)
        if url is None:
            return None
        url = "http://api.brain-map.org" + url
        interactive_download(url, cache_file)
        
    return cache_file

import json
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen
from .util import iter_download_with_resume


warehouse_url = 'http://iwarehouse'
warehouse_index_url = warehouse_url + '/api/v2/data/WellKnownFile/query.json?criteria=[path$il*synphys*]&num_rows=5000'

_index = None
def get_index():
    global _index
    if _index is None:
        req = urlopen(warehouse_index_url)
        data = json.load(req)
        _index = {rec['attachable_id']: rec['download_link'] for rec in data['msg']}
    return _index


def get_result_url(result_id):
    ind = get_index()
    return warehouse_url + ind.get(result_id)


def download_result_file(result_id, dest_file):
    url = get_result_url(result_id)
    assert url is not None
    for i,tot in iter_download_with_resume(url, dest_file, chunksize=100000000):
        print("%d / %d MB\r" % (i*1e-6, tot*1e-6))
    
    
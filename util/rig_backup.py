from __future__ import print_function
import os, sys
from multipatch_analysis.util import sync_file

args = sys.argv[1:]
if '--test' in args:
    test = True
    args.remove('--test')
else:
    test = False
    
source_path, dest_path = args
assert os.path.isdir(source_path)
assert os.path.isdir(dest_path)

for path, subdirs, files in os.walk(source_path):
    rel = os.path.relpath(path, source_path)
    for fname in files:
        src_file = os.path.join(path, fname)
        dst_file = os.path.join(dest_path, rel, fname)
        action = sync_file(src_file, dst_file, test=test)

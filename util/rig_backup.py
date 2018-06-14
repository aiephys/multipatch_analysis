from __future__ import print_function
import os, sys
from multipatch_analysis.util import sync_dir

args = sys.argv[1:]
if '--test' in args:
    test = True
    args.remove('--test')
else:
    test = False
    
source_path, dest_path = args
sync_dir(source_path, dest_path, test=test)

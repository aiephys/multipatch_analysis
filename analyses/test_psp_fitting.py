"""test psp fitting function using data from test directory
"""
import os
import numpy as np
from pprint import pprint
import json
import neuroanalysis.data
from multipatch_analysis.connection_detection import fit_psp

test_data_dir='test_psp_fit' # directory containing test data

test_data_files=[os.path.join(test_data_dir,f) for f in os.listdir(test_data_dir)] #list of test files
for file in test_data_files:
    test_dict=json.load(open(file)) # load test data
    avg_trace=neuroanalysis.data.Trace(data=np.array(test_dict['input']['data'])) # create Trace object
    psp_fits = fit_psp(avg_trace, 
                       sign=test_dict['input']['amp_sign'], 
                       yoffset=test_dict['input']['yoffset'], 
                       xoffset=test_dict['input']['xoffset'], 
                       amp=test_dict['input']['avg_amp'],
                       method=test_dict['input']['method'], 
                       stacked=test_dict['input']['stacked'], 
                       rise_time_mult_factor=test_dict['input']['rise_time_mult_factor'], 
                       fit_kws={'weight': test_dict['input']['weight']})                        

        
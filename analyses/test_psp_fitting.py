"""test psp fitting function using data from test directory
"""
import os
import numpy as np
from pprint import pprint
import json
import neuroanalysis.data

test_data_dir='test_psp_fit' # directory containing test data

test_data_files=[os.path.join(test_data_dir,f) for f in os.listdir(test_data_dir)] #list of test files
for file in test_data_files:
    test_dict=json.load(open(file)) # load test data
    avg_trace=neuroanalysis.data.Trace(data=np.array(test_dict['input']['data'])) # create Trace object
    psp_fits = fit_psp(avg_trace, 
                       sign=save_dict['input']['amp_sign'], 
                       yoffset=save_dict['input'][yoffset], 
                       xoffset=save_dict['input'][xoffset], 
                       amp=save_dict['input']['avg_amp'],
                       method=save_dict['input']['method'], 
                       stacked = save_dict['input']['stacled'], 
                       rise_time_mult_factor=10., 
                       fit_kws={'weights': weight})                        

        
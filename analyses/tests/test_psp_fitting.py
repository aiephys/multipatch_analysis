'''tests the psp_fits function.  If tests fail can use check_psp_fitting() to 
explore how fits are different from test data.  Often fitts are just off by
significant figures.  However, make sure all data is viewed as some fits can 
can change when others do not
'''
from multipatch_analysis.connection_detection import fit_psp
import os
import numpy as np
from pprint import pprint
import json
import neuroanalysis.data

test_data_dir = os.path.join(os.path.dirname(neuroanalysis.__file__), '..', 'test_data', 'test_psp_fit')                                                                                                                                      
    
def test_psp_fitting():
    """Test psp_fit function against data from test directory.  Note that this 
    test is highly sensitive.  If this test fails check_psp_fitting can be 
    used to investigate whether the differences are substantial. Many things
    can change the output of the fit slightly that would not be considered a real
    difference from a scientific perspective.  i.e. numbers off by a precision 
    of e-6.  One should look though the plots created by check_psp_fitting if there
    is any question.  Unexpected things such as the order of the parameters fed 
    to the function can create completely different fits.
    """
    plotting=True # specifies whether to make plots of fitting results
        
    test_data_files=[os.path.join(test_data_dir,f) for f in os.listdir(test_data_dir)] #list of test files
    for file in sorted(test_data_files):
#    for file in ['test_psp_fit/1492546902.92_2_6stacked.json']: order of parameters affects this fit
        print 'file', file
        test_dict=json.load(open(file)) # load test data
        avg_trace=neuroanalysis.data.Trace(data=np.array(test_dict['input']['data']), dt=test_dict['input']['dt']) # create Trace object
        psp_fits = fit_psp(avg_trace, 
                           xoffset=(14e-3, -float('inf'), float('inf')),
                           weight=np.array(test_dict['input']['weight']),
                           sign=test_dict['input']['amp_sign'], 
                           stacked=test_dict['input']['stacked'] 
                            )                        
        
        assert test_dict['out']['best_values']==psp_fits.best_values, \
            "Best values don't match. Run check_psp_fitting for more information"

        assert test_dict['out']['best_fit']==psp_fits.best_fit.tolist(), \
            "Best fit traces don't match. Run check_psp_fitting for more information"

        assert test_dict['out']['nrmse']==float(psp_fits.nrmse()), \
           "Nrmse doesn't match. Run check_psp_fitting for more information"

 

def check_psp_fitting():
    """Plots the results of the current fitting with the save fits and denotes 
    when there is a change. 
    """
    plotting=True # specifies whether to make plots of fitting results
        
    test_data_files=[os.path.join(test_data_dir,f) for f in os.listdir(test_data_dir)] #list of test files
    for file in sorted(test_data_files):
#    for file in ['test_psp_fit/1492546902.92_2_6stacked.json']: order of parameters affects this fit
        print 'file', file
        test_dict=json.load(open(file)) # load test data
        avg_trace=neuroanalysis.data.Trace(data=np.array(test_dict['input']['data']), dt=test_dict['input']['dt']) # create Trace object
        psp_fits = fit_psp(avg_trace, 
                           weight=np.array(test_dict['input']['weight']),
                           xoffset=(14e-3, -float('inf'), float('inf')),
                           sign=test_dict['input']['amp_sign'], 
                           stacked=test_dict['input']['stacked'] 
                            )                        
        
        change_flag=False
        if test_dict['out']['best_values']!=psp_fits.best_values:     
            print '  the best values dont match'
            print '\tsaved', test_dict['out']['best_values']
            print '\tobtained', psp_fits.best_values
            change_flag=True
            
        if test_dict['out']['best_fit']!=psp_fits.best_fit.tolist():
            print '  the best fit traces dont match'
            print '\tsaved', test_dict['out']['best_fit']
            print '\tobtained', psp_fits.best_fit.tolist()
            change_flag=True
        
        if test_dict['out']['nrmse']!=float(psp_fits.nrmse()):
            print '  the nrmse doesnt match'
            print '\tsaved', test_dict['out']['nrmse']
            print '\tobtained', float(psp_fits.nrmse())
            change_flag=True
            
        if plotting:
            import matplotlib.pylab as mplt
            fig=mplt.figure(figsize=(20,8))
            ax=fig.add_subplot(1,1,1)
            ax2=ax.twinx()
            ax.plot(avg_trace.time_values, psp_fits.data*1.e3, 'b', label='data')
            ax.plot(avg_trace.time_values, psp_fits.best_fit*1.e3, 'g', lw=5, label='current best fit')
            ax2.plot(avg_trace.time_values, test_dict['input']['weight'], 'r', label='weighting')
            if change_flag is True:
                ax.plot(avg_trace.time_values, np.array(test_dict['out']['best_fit'])*1.e3, 'k--', lw=5, label='original best fit')
                mplt.annotate('CHANGE', xy=(.5, .5), xycoords='figure fraction', fontsize=40)
            ax.legend()
            mplt.title(file + ', nrmse =' + str(psp_fits.nrmse()))
            mplt.show()
            
if __name__== "__main__":

#    check_psp_fitting() #use this to diagnose how fits differ from test data
    test_psp_fitting()

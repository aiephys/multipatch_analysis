from copy import deepcopy
import numpy as np
import scipy.signal
import pyqtgraph as pg

from .data import MultiPatchProbe, Analyzer, PulseStimAnalyzer
from . import qc
from neuroanalysis.stats import ragged_mean
from neuroanalysis.data import Trace, TraceList
from neuroanalysis.fitting import StackedPsp, Psp
from neuroanalysis.ui.plot_grid import PlotGrid
from neuroanalysis.filter import bessel_filter

class BaselineDistributor(Analyzer):
    """Used to find baseline regions in a trace and distribute them on request.
    """
    def __init__(self, rec):
        self._attach(rec)
        self.rec = rec
        self.baselines = rec.baseline_regions
        self.ptr = 0

    def get_baseline_chunk(self, duration=20e-3):
        """Return the (start, stop) indices of a chunk of unused baseline with the
        given duration.
        """
        dt = self.rec['primary'].dt
        chunk_size = int(np.round(duration / dt))
        while True:
            if len(self.baselines) == 0:
                return None
            start, stop = self.baselines[0]
            chunk_start = max(start, self.ptr)
            chunk_stop = chunk_start + chunk_size
            if chunk_stop <= stop:
                self.ptr = chunk_stop
                return chunk_start, chunk_stop
            else:
                self.baselines.pop(0)


class MultiPatchSyncRecAnalyzer(Analyzer):
    """Used for analyzing two or more synchronous patch clamp recordings where
    spikes are evoked in at least one and synaptic responses are recorded in
    others.
    """
    def __init__(self, srec):
        self._attach(srec)
        self.srec = srec

    def get_spike_responses(self, pre_rec, post_rec, align_to='pulse', pre_pad=10e-3, require_spike=True):
        """Given a pre- and a postsynaptic recording, return a structure
        containing evoked responses.
        
            [{pulse_n, pulse_ind, spike, response, baseline}, ...]
        
        """
        # detect presynaptic spikes
        pulse_stim = PulseStimAnalyzer.get(pre_rec)
        spikes = pulse_stim.evoked_spikes()
        
        if not isinstance(pre_rec, MultiPatchProbe):
            # this does not look like the correct kind of stimulus; bail out
            # Ideally we can make this agnostic to the exact stim type in the future,
            # but for now we rely on the delay period between pulses 8 and 9 to get
            # a baseline measurement.
            return []

        dt = post_rec['primary'].dt
        
        # Select ranges to extract from postsynaptic recording
        result = []
        for i,pulse in enumerate(spikes):
            pulse = pulse.copy()
            spike = pulse['spike']
            if require_spike and spike is None:
                continue
            
            if align_to == 'spike':
                # start recording window at the rising phase of the presynaptic spike
                pulse['rec_start'] = spike['rise_index'] - int(pre_pad / dt)
            elif align_to == 'pulse':
                # align to pulse onset
                pulse['rec_start'] = pulse['pulse_ind'] - int(pre_pad / dt)
            pulse['rec_start'] = max(0, pulse['rec_start'])
            
            # get times of nearby pulses
            prev_pulse = None if i == 0 else spikes[i-1]['pulse_ind'] + spikes[i-1]['pulse_len']
            this_pulse = pulse['pulse_ind']
            next_pulse = None if i+1 >= len(spikes) else spikes[i+1]['pulse_ind']

            max_stop = pulse['rec_start'] + int(50e-3 / dt)
            if next_pulse is not None:
                # truncate window early if there is another pulse
                pulse['rec_stop'] = min(max_stop, next_pulse)
            else:
                # otherwise, stop 50 ms later
                pulse['rec_stop'] = max_stop
            
            # Extract data from postsynaptic recording
            pulse['response'] = post_rec['primary'][pulse['rec_start']:pulse['rec_stop']]
            assert len(pulse['response']) > 0

            # Extract presynaptic spike and stimulus command
            pulse['pre_rec'] = pre_rec['primary'][pulse['rec_start']:pulse['rec_stop']]
            pulse['command'] = pre_rec['command'][pulse['rec_start']:pulse['rec_stop']]

            # select baseline region between 8th and 9th pulses
            baseline_dur = int(100e-3 / dt)
            stop = spikes[8]['pulse_ind']
            start = stop - baseline_dur
            pulse['baseline'] = post_rec['primary'][start:stop]
            pulse['baseline_start'] = start
            pulse['baseline_stop'] = stop

            # Add minimal QC metrics for excitatory and inhibitory measurements
            pulse_window = [pulse['rec_start'], pulse['rec_stop']]
            n_spikes = 0 if spike is None else 1  # eventually should check for multiple spikes
            adj_pulse_times = []
            if prev_pulse is not None:
                adj_pulse_times.append((prev_pulse - this_pulse) * dt)
            if next_pulse is not None:
                adj_pulse_times.append((next_pulse - this_pulse) * dt)
            pulse['ex_qc_pass'], pulse['in_qc_pass'] = qc.pulse_response_qc_pass(post_rec=post_rec, window=pulse_window, n_spikes=n_spikes, adjacent_pulses=adj_pulse_times)

            assert len(pulse['baseline']) > 0

            result.append(pulse)
        
        return result

    def get_pulse_response(self, pre_rec, post_rec, first_pulse=0, last_pulse=-1):
        pulse_stim = PulseStimAnalyzer.get(pre_rec)
        spikes = pulse_stim.evoked_spikes()
        
        if len(spikes) < 10:
            return None
        
        dt = post_rec['primary'].dt
        start = spikes[first_pulse]['pulse_ind'] - int(20e-3 / dt)
        stop = spikes[last_pulse]['pulse_ind'] + int(50e-3 / dt)
        return post_rec['primary'][start:stop]

    def stim_params(self, pre_rec):
        return PulseStimAnalyzer.get(pre_rec).stim_params()
        
    def get_train_response(self, pre_rec, post_rec, start_pulse, stop_pulse, padding=(-10e-3, 50e-3)):
        """Return the part of the post-synaptic recording during a range of pulses,
        along with a baseline chunk
        """
        pulse_stim = PulseStimAnalyzer.get(pre_rec)
        pulses = [p[0] for p in pulse_stim.pulses() if p[2] > 0]
        
        post_trace = post_rec['primary']
        pre_trace = pre_rec['primary']
        stim_command = pre_rec['command']

        dt = post_trace.dt
        start = pulses[start_pulse] + int(padding[0]/dt)
        stop = pulses[stop_pulse] + int(padding[1]/dt)
        assert start > 0
        
        response = post_trace[start:stop]
        pre_spike = pre_trace[start:stop]
        command = stim_command[start:stop]
        
        bstop = pulses[8] - int(10e-3/dt)
        bstart = bstop - int(200e-3/dt)
        baseline = post_trace[bstart:bstop]
        
        return response, baseline, pre_spike, command

    def find_artifacts(self, rec, freq=None, pos_threshold=-10e-3, neg_threshold=-100e-3):
        """Find contaminating artifacts in recording such as spontaneous spikes or electrical noise, return True
        state if anything is found"""

        artifacts = False
        # find anything that crosses a max threshold of -10mV
        if np.any(rec >= pos_threshold):
            artifacts = True
        if np.any(rec <= neg_threshold):
            artifacts = True
        return artifacts
        


class MultiPatchExperimentAnalyzer(Analyzer):
    """
    loop over all sweeps (presynaptic)
        ignore sweeps with high induction frequency
        detect presynaptic spikes
        loop over all other sweeps in the same recording (postsynaptic)
            ignore sweeps that fail QC
            collect data following each stimulus
            collect baseline noise data
                use noise.std() / sqrt(n) to estimate noise for averaged traces
                use mode(noise) to correct for baseline
            collect mode / holding
    Average together all collected sweep chunks from the same cell, and of the same holding and clamp mode
        all pulses
        pulses 1, 2, 9, 10
        pulses 7, 8, 11, 12
    Fit result to psp
    Compare fit amplitude to noise floor
    Z-score fit NRMSE against fits from pre-stimulus data
        correct for multiple comparisons?
    Additional metric to detect low release probability connections?
    """
    def __init__(self, expt):
        self._attach(expt)
        self.expt = expt
        self._all_spikes = None

    def get_evoked_responses(self, pre_id, post_id, clamp_mode='ic', stim_filter=None, min_duration=None, pulse_ids=None):
        """Return all evoked responses from device pre_id to post_id with the given
        clamp mode and stimulus conditions.
        
        All traces are *downsampled* to the minimum sample rate across the set
        of returned responses.
        
        Returns a list of (response, baseline) pairs. 
        """
        all_spikes = self._all_evoked_responses()
        responses = EvokedResponseGroup(pre_id, post_id)
        if pre_id not in all_spikes or post_id not in all_spikes[pre_id]:
            return responses
        
        for rec in all_spikes[pre_id][post_id]:
            
            # do filtering here:
            pre_rec = rec['pre_rec']
            post_rec = rec['post_rec']
            
            if post_rec.clamp_mode != clamp_mode:
                continue
            
            if stim_filter is not None:
                stim_name = pre_rec.stimulus.description
                if stim_filter not in stim_name:
                    continue
            
            for spike in rec['spikes']:
                if spike['spike'] is None:
                    continue
                if pulse_ids is not None and spike['pulse_n'] not in pulse_ids:
                    continue
                resp = spike['response']
                if resp.duration >= min_duration:
                    responses.add(resp.copy(t0=0), spike['baseline'])
        
        return responses
 
    def get_evoked_response_matrix(self, **kwds):
        """Returned evoked responses for all pre/post pairs
        """
        devs = self.list_devs()

        all_responses = {}
        rows = set()
        cols = set()
        for i, dev1 in enumerate(devs):
            for j, dev2 in enumerate(devs):
                if dev1 == dev2:
                    continue
                resp = self.get_evoked_responses(dev1, dev2, **kwds)
                if len(resp) > 0:
                    rows.add(dev1)
                    cols.add(dev2)
                all_responses[(dev1, dev2)] = resp
        rows = sorted(list(rows))
        cols = sorted(list(cols))

        return all_responses, rows, cols

    def _all_evoked_responses(self):
        
        # loop over all sweeps (presynaptic)
        if self._all_spikes is None:
            all_spikes = {}
            for srec in self.expt.contents:
                mp_analyzer = MultiPatchSyncRecAnalyzer.get(srec)
                
                for pre_rec in srec.recordings:
                    if not isinstance(pre_rec, MultiPatchProbe):
                        continue
                    pre_id = pre_rec.device_id
                    all_spikes.setdefault(pre_id, {})
                    # todo: ignore sweeps with high induction frequency
                    
                    for post_rec in srec.recordings:
                        if post_rec is pre_rec:
                            continue
                        post_id = post_rec.device_id
                        all_spikes[pre_id].setdefault(post_id, [])
                        spikes = {
                            'spikes': mp_analyzer.get_spike_responses(pre_rec, post_rec),
                            'pre_rec': pre_rec,
                            'post_rec': post_rec,
                        }
                        all_spikes[pre_id][post_id].append(spikes)
            self._all_spikes = all_spikes
        return self._all_spikes

    def list_devs(self):
        return self._all_evoked_responses().keys()



class NDDict(object):
    """N-dimensional dictionary.
    """
    def __init__(self, ndim):
        self._ndim = int(ndim)
        self._data = {}
        self._keys = [set() for i in range(ndim)]

    @property
    def ndim(self):
        return self._ndim

    def __setitem__(self, item, value):
        ndim = self.ndim
        assert len(item) == ndim
        data = self._data
        for i,key in enumerate(item[:-1]):
            data = data.setdefault(key, {})
            self.keys()[i].add(key)

        data[item[-1]] = value
        self.keys()[ndim-1].add(item[-1])

    def keys(self):
        return self._keys



def plot_response_averages(expt, show_baseline=False, **kwds):
    analyzer = MultiPatchExperimentAnalyzer.get(expt)
    devs = analyzer.list_devs()

    # First get average evoked responses for all pre/post pairs
    responses, rows, cols = analyzer.get_evoked_response_matrix(**kwds)

    # resize plot grid accordingly
    plots = PlotGrid()
    plots.set_shape(len(rows), len(cols))
    plots.show() 
    
    ranges = [([], []), ([], [])]
    points = []

    # Plot each matrix element with PSP fit
    for i, dev1 in enumerate(rows):
        for j, dev2 in enumerate(cols):
            # select plot and hide axes
            plt = plots[i, j]
            if i < len(devs) - 1:
                plt.getAxis('bottom').setVisible(False)
            if j > 0:
                plt.getAxis('left').setVisible(False)

            if dev1 == dev2:
                plt.getAxis('bottom').setVisible(False)
                plt.getAxis('left').setVisible(False)
                continue
            
            # adjust axes / labels
            plt.setXLink(plots[0, 0])
            plt.setYLink(plots[0, 0])
            plt.addLine(x=10e-3, pen=0.3)
            plt.addLine(y=0, pen=0.3)
            plt.setLabels(bottom=(str(dev2), 's'))
            if kwds.get('clamp_mode', 'ic') == 'ic':
                plt.setLabels(left=('%s' % dev1, 'V'))
            else:
                plt.setLabels(left=('%s' % dev1, 'A'))

            
            # print "==========", dev1, dev2
            avg_response = responses[(dev1, dev2)].bsub_mean()
            if avg_response is not None:
                avg_response.t0 = 0
                t = avg_response.time_values
                y = bessel_filter(Trace(avg_response.data, dt=avg_response.dt), 2e3).data
                plt.plot(t, y, antialias=True)

                # fit!                
                #fit = responses[(dev1, dev2)].fit_psp(yoffset=0, mask_stim_artifact=(abs(dev1-dev2) < 3))
                #lsnr = np.log(fit.snr)
                #lerr = np.log(fit.nrmse())
                
                #color = (
                    #np.clip(255 * (-lerr/3.), 0, 255),
                    #np.clip(50 * lsnr, 0, 255),
                    #np.clip(255 * (1+lerr/3.), 0, 255)
                #)

                #plt.plot(t, fit.best_fit, pen=color)
                ## plt.plot(t, fit.init_fit, pen='y')

                #points.append({'x': lerr, 'y': lsnr, 'brush': color})

                #if show_baseline:
                    ## plot baseline for reference
                    #bl = avg_response.meta['baseline'] - avg_response.meta['baseline_med']
                    #plt.plot(np.arange(len(bl)) * avg_response.dt, bl, pen=(0, 100, 0), antialias=True)

                # keep track of data range across all plots
                ranges[0][0].append(y.min())
                ranges[0][1].append(y.max())
                ranges[1][0].append(t[0])
                ranges[1][1].append(t[-1])

    plots[0,0].setYRange(min(ranges[0][0]), max(ranges[0][1]))
    plots[0,0].setXRange(min(ranges[1][0]), max(ranges[1][1]))

    # scatter plot of SNR vs NRMSE
    plt = pg.plot()
    plt.setLabels(left='ln(SNR)', bottom='ln(NRMSE)')
    plt.plot([p['x'] for p in points], [p['y'] for p in points], pen=None, symbol='o', symbolBrush=[pg.mkBrush(p['brush']) for p in points])
    # show threshold line
    line = pg.InfiniteLine(pos=[0, 6], angle=180/np.pi * np.arctan(1))
    plt.addItem(line, ignoreBounds=True)

    return plots


class EvokedResponseGroup(object):
    """A group of similar synaptic responses.

    This is intended to be used as a container for many repeated responses evoked from
    a single pre/postsynaptic pair. It provides methods for computing the average,
    baseline-subtracted response and for fitting the average to a curve.
    """
    def __init__(self, pre_id=None, post_id=None, **kwds):
        self.pre_id = pre_id
        self.post_id = post_id
        self.kwds = kwds
        self.responses = []
        self.baselines = []
        self.spikes = []
        self.commands = []
        self._bsub_mean = None

    def add(self, response, baseline, pre_spike=None, stim_command=None):
        self.responses.append(response)
        self.baselines.append(baseline)
        self.spikes.append(pre_spike)
        self.commands.append(stim_command)
        self._bsub_mean = None

    def __len__(self):
        return len(self.responses)

    def bsub_mean(self):
        """Return a baseline-subtracted, average evoked response trace between two cells.

        All traces are downsampled to the minimum sample rate in the set.
        """
        if len(self) == 0:
            return None

        if self._bsub_mean is None:
            responses = self.responses
            baselines = self.baselines
            
            # downsample all traces to the same rate
            # yarg: how does this change SNR?

            avg = TraceList([r.copy(t0=0) for r in responses]).mean()
            avg_baseline = TraceList([b.copy(t0=0) for b in baselines]).mean().data

            # subtract baseline
            baseline = np.median(avg_baseline)
            bsub = avg.data - baseline
            result = avg.copy(data=bsub)
            assert len(result.time_values) == len(result)

            # Attach some extra metadata to the result:
            result.meta['baseline'] = avg_baseline
            result.meta['baseline_med'] = baseline
            if len(avg_baseline) == 0:
                result.meta['baseline_std'] = None
            else:
                result.meta['baseline_std'] = scipy.signal.detrend(avg_baseline).std()

            self._bsub_mean = result

        return self._bsub_mean

    def mean(self):
        if len(self) == 0:
            return None
        return TraceList(self.responses).mean()

    def fit_psp(self, **kwds):
        response = self.bsub_mean()
        if response is None:
            return None
        return fit_psp(response, **kwds)
    
def test_create_all_fit_param_combos():
    """Tests for changes in create_all_fit_param_combos module.
    """
    input={'decay_tau': (0.05, 0.005, 0.5), 
           'rise_time': (0.005, 0.0005, 0.05), 
           'yoffset': (0, -float('inf'), float('inf')), 
           'rise_power': ([1, 2], 'fixed'), 
           'amp': (0.0002, 0, 0.1), 
           'xoffset': ([0.011, 0.014, 0.015], -float('inf'), float('inf'))}
    output=[{'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (1, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.011, -float('inf'), float('inf'))}, 
            {'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (2, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.011, -float('inf'), float('inf'))}, 
            {'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (1, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.014, -float('inf'), float('inf'))}, 
            {'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (2, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.014, -float('inf'), float('inf'))}, 
            {'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (1, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.015, -float('inf'), float('inf'))}, 
            {'decay_tau': (0.05, 0.005, 0.5), 'rise_time': (0.005, 0.0005, 0.05), 'yoffset': (0, -float('inf'), float('inf')), 'rise_power': (2, 'fixed'), 'amp': (0.0002, 0, 0.1), 'xoffset': (0.015, -float('inf'), float('inf'))}]

    test_out=create_all_fit_param_combos(input)
    
    assert test_out == output

def create_all_fit_param_combos(base_params):
    '''Convert the parameters fed to fit_psp into a list of all possible parameter 
    dictionaries to be fed to PSP() or stackedPSP() for fitting. 
    
    Parameters 
    ----------
    base_params: dictionary
        Each value in the key:value dictionary pair must be a tuple.
        In general the structure of the tuple is of the form, 
        (initial conditions, lower boundary, higher boundary).
        The initial conditions can be either a number or a list 
        of numbers specifying several initial conditions.  The 
        initial condition may also be fixed by replacing the lower 
        higher boundary combination with 'fixed'.
        Note that technically the value of a key:value pair could 
        be a single string (this ends up being evaluated as an 
        expression by lmfit later on).
    
    Returns
    -------
    param_dict_list: list of dictionaries
        Each dictionary contains parameter inputs for one fit_psp run.
        
    Examples:    
    base_params[amplitude]=(10, 0, 20)
    base_params[amplitude]=(10, 'fixed')
    base_params[amplitude]=([5,10, 20], 0, 20)
    base_params[amplitude]=([5,10, 20], 'fixed')
    
    '''
    # need to create all combinations of the initial conditions
    param_dict_list = [{}] #initialize list
    for key, value in base_params.iteritems():
        if isinstance(value[0], list):
            temp_list=[]
            for init_cond in value[0]: #for each initial condition
                temp=[pdl.copy() for pdl in param_dict_list] #copy each parameter dictionary in the list so they do not point back to the original dictionary
                for t in temp:  #for each dictionary in the list 
                    t[key]=tuple([init_cond] +list(value[1:])) #add the key and single initial condition pair
                temp_list=temp_list+temp
            param_dict_list=list(temp_list) #Note this works because the dictionaries are already made immutable above
        else: 
            for pdl in param_dict_list:
                pdl[key]=value
    
    return param_dict_list

def fit_psp(response, 
            mode='ic', 
            sign='any', #Note this will not be used if *amp* input is specified
            method='leastsq', 
            fit_kws=None, 
            stacked=True,
            rise_time_mult_factor=10., #Note this will not be used if *rise_time* input is specified 
            weight='default',
            amp_ratio='default', 
            # the following are parameters that can be fit 
                amp='default',
                decay_tau='default',
                rise_power='default',
                rise_time='default',
                exp_amp='default',
                xoffset='default', 
                yoffset='default',
            ):
    """Fit psp. function to the equation 
    
    This function make assumptions about where the cross talk happens as traces 
    have been aligned to the pulses and response.t0 has been set to zero 
    
    Parameters
    ----------
    response : neuroanalysis.data.Trace class
        Contains data on trace waveform.
    mode : string
        either 'ic' for current clamp or 'vc' for voltage clamp
    sign : string
        Specifies the sign of the PSP deflection.  Must be '+', '-', or any. If *amp* 
        is specified, value will be irrelevant.
    method : string 
        Method lmfit uses for optimization
    rise_time_mult_factor: float
        Parameter that goes into the default calculation rise time.  
        Note that if an input for *rise_time* is provided this input
        will be irrelevant.
    stacked : True or False
        If true, use the StackedPsp function which assumes there is still left
        over voltage decay from previous events.  If False, use Psp function
        which assumes the region of the waveform before the event is at baseline.
    fit_kws : dictionary
        Additional key words that are fed to lmfit
    exp_amp : string
        function that is fed to lmfit
    The parameters below are fed to the psp function. Each value in the 
        key:value dictionary pair must be a tuple.
        In general the structure of the tuple is of the form, 
        (initial conditions, lower boundary, higher boundary).
        The initial conditions can be either a number or a list 
        of numbers specifying several initial conditions.  The 
        initial condition may also be fixed by replacing the lower 
        higher boundary combination with 'fixed'.    
        Examples:    
            amplitude=(10, 0, 20)
            amplitude=(10, 'fixed')
            amplitude=([5,10, 20], 0, 20)
            amplitude=([5,10, 20], 'fixed') 
        xoffset : scalar
            Horizontal shift between begin (positive shifts to the right)
        yoffset : scalar
            Vertical offset
        rise_time : scalar
            Time from beginning of psp until peak
        decay_tau : scalar
            Decay time constant
        amp : scalar
            The peak value of the psp
        rise_power : scalar
            Exponent for the rising phase; larger values result in a slower activation 
        amp_ratio : scalar 
            if *stacked* this is used to set up the ratio between the 
            residual decay amplitude and the height of the PSP.
    
    Returns
    -------
    fit: lmfit.model.ModelResult
        Best fit
    """           
    
    # extracting these for ease of use
    t=response.time_values
    y=response.data
    dt = response.dt
    
    # set initial conditions depending on whether in voltage or current clamp
    # note that sign of these will automatically be set later on based on the 
    # the *sign* input
    if mode == 'ic':
        amp_init = .2e-3
        amp_max = 100e-3
        rise_time_init = 5e-3
        decay_tau_init = 50e-3
    elif mode == 'vc':
        amp_init = 20e-12
        amp_max = 500e-12
        rise_time_init = 1e-3
        decay_tau_init = 4e-3
    else:
        raise ValueError('mode must be "ic" or "vc"')

    # Set up amplitude initial values and boundaries depending on whether *sign* are positive or negative
    if sign == '-':
        amps = (-amp_init, -amp_max, 0)
    elif sign == '+':
        amps = (amp_init, 0, amp_max)
    elif sign =='any':
        warnings.warn("You are not specifying the predicted sign of your psp.  This may slow down or mess up fitting")
        amps = (0, -amp_max, amp_max)
    else:
        raise ValueError('sign must be "+", "-", or "any"')
        
    # initial condition, lower boundry, upper boundry    
    base_params = {
        'xoffset': (14e-3, -float('inf'), float('inf')),
        'yoffset': (0, -float('inf'), float('inf')),
        'rise_time': (rise_time_init, rise_time_init/rise_time_mult_factor, rise_time_init*rise_time_mult_factor),
        'decay_tau': (decay_tau_init, decay_tau_init/10., decay_tau_init*10.),
        'rise_power': (2, 'fixed'),
        'amp': amps
    }
    
    # specify fitting function and set up conditions
    if not isinstance(stacked, bool):
        raise Exception("Stacked must be True or False")
    if stacked:
        psp = StackedPsp()
        base_params.update({
            #TODO: figure out the bounds on these
            'exp_amp': 'amp * amp_ratio',
            'amp_ratio': (0, -100, 100),
        })  
    else:
        psp = Psp()
    
    # override defaults with input
    for bp in base_params.keys():
        if eval(bp)!='default':
            base_params[bp]=eval(bp)
    
    # set weighting that 
    if weight == 'default': #use default weighting
        # THIS CODE IS DEPENDENT ON THE DATA BEING INPUT IN A CERTAIN WAY THAT IS NOT TESTED
        weight = np.ones(len(y))*10.  #set everything to ten initially
        weight[int(10e-3/dt):int(12e-3/dt)] = 0.   #area around stim artifact
        weight[int(12e-3/dt):int(19e-3/dt)] = 30.  #area around steep PSP rise 
    elif weight is False: #do not weight any part of the stimulus
        weight = np.ones(len(y))
    elif 'weight' in vars():  #works if there is a value specified in weight
        if len(weight) != len(y):
            raise Exception('the weight and array vectors are not the same length') 
    
    # arguement to be passed through to fitting function
    fit_kws={'weights': weight}   

    # convert initial parameters into a list of dictionaries to be consumed by psp.fit()        
    param_dict_list= create_all_fit_param_combos(base_params)

    # cycle though different parameters sets and chose best one
    best_fit = None
    best_score = None
    for p in param_dict_list:
        fit = psp.fit(y, x=t, params=p, fit_kws=fit_kws, method=method)
        err = np.sum(fit.residual**2)  # note: using this because normalized (nrmse) is not necessary to comparing fits within the same data set
        if best_fit is None or err < best_score:
            best_fit = fit
            best_score = err
    fit = best_fit

    # nrmse = fit.nrmse()
    if 'baseline_std' in response.meta:
        fit.snr = abs(fit.best_values['amp']) / response.meta['baseline_std']
        fit.err = fit.nrmse() / response.meta['baseline_std']

    return fit


def detect_connections(expt):
    analyzer = MultiPatchExperimentAnalyzer.get(expt)

    # First get average evoked responses for all pre/post pairs with long decay time
    all_responses, rows, cols = analyzer.get_evoked_response_matrix(clamp_mode='ic', min_duration=16e-3)

    for pre_id in rows:
        for post_id in cols:
            try:
                response = all_responses[(pre_id, post_id)]
                if len(response) == 0:
                    continue
            except KeyError:
                continue

            # fit average to extract PSP decay
            fit = response.fit_psp(yoffset=0)

            # make connectivity call
            lsnr = np.log(fit.snr)
            lnrmse = np.log(fit.nrmse())
            if lsnr > lnrmse + 6:
                print "Connection:", pre_id, post_id, fit.snr, fit.nrmse()



from datetime import datetime
import database as db
from neuroanalysis.baseline import float_mode
from neuroanalysis.data import PatchClampRecording
from lims import specimen_info
from allensdk_internal.core import lims_utilities as lims
from data import MultipatchExperiment
from connection_detection import PulseStimAnalyzer, MultiPatchSyncRecAnalyzer
import config



class SliceSubmission(object):
    """Used to submit a new slice entry to the synphys DB.
    """
    def __init__(self, dh):
        self.dh = dh
        
        self._fields = None

    @property
    def fields(self):
        if self._fields is None:
            info = self.dh.info()
            
            # pull some metadata from LIMS
            sid = info['specimen_ID']
            limsdata = specimen_info(sid)

            self._fields = {
                'acq_timestamp': datetime.fromtimestamp(info['__timestamp__']),
                'species': limsdata['organism'],
                'age': limsdata['age'],
                'genotype': limsdata['genotype'],
                'orientation': limsdata['plane_of_section'],
                'surface': limsdata['exposed_surface'],
                'hemisphere': limsdata['hemisphere'],
                'quality': info.get('slice quality'),
                'slice_time': info.get('slice time'),
                'slice_conditions': {},
                'lims_specimen_name': sid,
                'original_path': '%s:%s' % (config.rig_name, self.dh.name()),
                'submission_data': None,
            }
        return self._fields

    def check(self):
        warnings = []
        errors = []
        fields = self.fields
        
        # TODO: Add a lot more checking here..
        
        if fields['surface'] not in ['right', 'left']:
            warnings.append("Warning: slice surface '%s' should have been 'right' or 'left'" % surface)
        
        return errors, warnings
        
    def create(self):
        if len(self.check()[0]) > 0:
            raise Exception("Submission has errors; see SliceSubmission.check()")
        data = self.fields
        sl = db.Slice(**data)
        return sl
        
    def submit(self):
        session = db.Session()
        sl = self.create()
        session.add(sl)
        session.commit()
        

def submit_slice(data):
    """Submit information about a new slice to the internal analysis DB.
    
        data = {
            'specimen_id': <LIMS specimen ID>,
            'original_path': <original file location eg: \\RIG\\D\...\slice_000>,
            'acquisition_uid': <unique ID chosen by acquisition system>,
            'surface': <recorded surface (medial/lateral)>,
            'image_files': {
                'slice anatomy': [image file names],
                'slice quality': [image file names],
            }
        }
    """
    
    
class ExperimentSubmission(object):
    """Used to submit a new experiment entry to the synphys DB.
    
    This causes several tables to be populated: experiment, sync_rec, recording,
    stim_pulse, stim_spike, pulse_response, baseline
    """
    def __init__(self, dh, nwb_file):
        self.dh = dh
        self.nwb_file = nwb_file
        self._fields = None

    @property
    def fields(self):
        if self._fields is None:
            info = self.dh.info()
            
            slice_dir = self.dh.parent()
            ts = datetime.fromtimestamp(slice_dir.info()['__timestamp__'])
            slice_id = db.slice_from_timestamp(ts).id
            
            expt_dir = slice_dir.parent()
            expt_info = expt_dir.info()

            temp = expt_info.get('temperature')
            if temp is not None:
                temp = float(temp.rstrip(' C'))

            self._fields = {
                'original_path': '%s:%s' % (config.rig_name, self.dh.name()),
                'acq_timestamp': datetime.fromtimestamp(info['__timestamp__']),
                'slice_id': slice_id,
                'target_region': expt_info.get('region'),
                'internal': expt_info.get('internal'),
                'acsf': expt_info.get('acsf'),
                'target_temperature': temp,
                
            }
        return self._fields

    def check(self):
        warnings = []
        errors = []
        fields = self.fields
        
        # TODO: Add a lot more checking here..
        
        return errors, warnings
        
    def create(self, session):
        if len(self.check()[0]) > 0:
            raise Exception("Submission has errors; see SiteSubmission.check()")
        
        # Create entry in experiment table
        data = self.fields
        expt = db.Experiment(**data)
        session.add(expt)
        
        # Load NWB file and create data entries
        nwb = MultipatchExperiment(self.nwb_file.name())

        for srec in nwb.contents:
            srec_entry = db.SyncRec(sync_rec_key=srec.key, experiment=expt)
            session.add(srec_entry)
            
            rec_entries = {}
            all_pulse_entries = {}
            for rec in srec.recordings:
                psa = PulseStimAnalyzer.get(rec)
                ind_freq, recovery_delay = psa.stim_params()
                
                # import all recordings
                rec_entry = db.Recording(
                    sync_rec=srec_entry,
                    device_key=rec.device_id, 
                    start_time=rec.start_time,
                )
                session.add(rec_entry)
                rec_entries[rec.device_id] = rec_entry
                
                # import patch clamp recording information
                if isinstance(rec, PatchClampRecording):
                    pcrec_entry = db.PatchClampRecording(
                        recording=rec_entry,
                        clamp_mode=rec.clamp_mode,
                        patch_mode=rec.patch_mode,
                        stim_name=rec.meta['stim_name'],
                        baseline_potential=rec.baseline_potential,
                        baseline_current=rec.baseline_current,
                        baseline_rms_noise=rec.baseline_rms_noise,
                    )
                    session.add(pcrec_entry)

                # import test pulse information
                tp = rec.nearest_test_pulse
                if tp is not None:
                    tp_entry = db.TestPulse(
                        patch_clamp_recording=pcrec_entry,
                        start_index=tp.indices[0],
                        stop_index=tp.indices[1],
                        baseline_current=tp.baseline_current,
                        baseline_potential=tp.baseline_potential,
                        access_resistance=tp.access_resistance,
                        input_resistance=tp.input_resistance,
                        capacitance=tp.capacitance,
                        time_constant=tp.time_constant,
                    )
                    session.add(tp_entry)
                
                # import presynaptic stim pulses
                pulses = psa.pulses()
                
                pulse_entries = {}
                all_pulse_entries[rec.device_id] = pulse_entries
                
                for i,pulse in enumerate(pulses):
                    if i == 0 and rec.has_inserted_test_pulse:
                        continue

                    pulse_entry = db.StimPulse(
                        recording=rec_entry,
                        pulse_number=i,
                        onset_index=pulse[0],
                        amplitude=pulse[2],
                        length=pulse[1]-pulse[0],
                    )
                    session.add(pulse_entry)
                    pulse_entries[i] = pulse_entry

                # import presynaptic evoked spikes
                spikes = psa.evoked_spikes()
                for i,sp in enumerate(spikes):
                    pulse = pulse_entries[sp['pulse_n']]
                    if sp['spike'] is not None:
                        extra = sp['spike']
                        pulse.n_spikes = 1
                    else:
                        extra = {}
                        pulse.n_spikes = 0
                    
                    spike_entry = db.StimSpike(
                        recording=rec_entry,
                        pulse=pulse,
                        **extra
                    )
                    session.add(spike_entry)
                
            # import postsynaptic responses
            mpa = MultiPatchSyncRecAnalyzer(srec)
            for pre_dev in srec.devices:
                for post_dev in srec.devices:
                    # get all responses, regardless of the presence of a spike
                    responses = mpa.get_spike_responses(srec[pre_dev], srec[post_dev], align_to='pulse', require_spike=False)
                    for resp in responses:
                        base_entry = db.Baseline(
                            recording=rec_entries[post_dev],
                            start_index=resp['baseline_start'],
                            stop_index=resp['baseline_stop'],
                            data=resp['baseline'].downsample(f=50000).data,
                            mode=float_mode(resp['baseline'].data),
                        )
                        session.add(base_entry)
                        resp_entry = db.PulseResponse(
                            recording=rec_entries[post_dev],
                            stim_pulse=all_pulse_entries[pre_dev][resp['pulse_n']],
                            baseline=base_entry,
                            start_index=resp['rec_start'],
                            stop_index=resp['rec_stop'],
                            data=resp['response'].downsample(f=50000).data,
                        )
                        session.add(resp_entry)
        return expt
        
    def submit(self):
        session = db.Session()
        exp = self.create(session)
        session.commit()



def submit_experiment(data):
    """Submit a new experiment to the internal analysis DB.
    
        data = {
            'slice_specimen_id': <LIMS specimen ID of parent slice>,
            'original_site_path': <original file location eg: \\RIG\\D\...\site_000>,
            'acquisition_uid': <unique ID chosen by acquisition system>,
            'nwb_file': <subpath to NWB file>,
            'images': {
                'recording site': [image file names],
            }
            'electrodes': 
        }
    
    
    Causes new structures to be generated in DB:
    
    * an entry in experiment table
    * entries in electrode, cell, pair tables
    * entries in syncrec, recording, trace tables 
      (these just enumerate the structures in the NWB)
    * entries in stim_pulse, stim_spike, pulse_response, and baseline tables
      (these are extracted from NWB)
    
    """
    
    prof = Profiler(disabled=True, delayed=False)
    
    if expt.expt_id in expt_index:
        print("Cached: %s" % expt)
        raise NotImplementedError()
    
    prof.mark('start')
    expt_id = len(expt_rows)
    expt_rows.append({'id': expt_id, 'expt_key': expt.expt_id, 'internal_id': -1,
        'acsf_id': -1, 'temperature': np.nan, 'age': expt.age, 'genotype': None,
        'date': expt.date})
    expt_index[expt.expt_id] = expt_id
    
    cell_ids = {}
    for cell in expt.cells.values():
        cell_id = len(cell_rows)
        # mapping from experiment's internal ID for this cell to global cell ID 
        cell_ids[cell.cell_id] = cell_id
        cell_rows.append({'id': cell_id, 'expt_id': expt_id, 
            'device_key': cell.cell_id, 'cre_type': cell.cre_type,
            'pass_qc': cell.pass_qc, 'position': cell.position,
            'depth': cell.depth})
    prof.mark('cells')


    expt_data = expt.data
    for srec in expt_data.contents:
        srec_id = len(srec_rows)
        srec_rows.append({'id': srec_id, 'expt_id': expt_id,
            'sync_rec_key': srec.key})
        rec_key_id_map = {}
        pulse_key_n_id_map = {}
        for rec in srec.recordings:
            rec_id = len(rec_rows)
            rec_key_id_map[rec.device_id] = rec_id
            tp_id = len(tp_rows)
            cell_id = cell_ids[rec.device_id + 1]
            psa = PulseStimAnalyzer.get(rec)
            ind_freq, recovery_delay = psa.stim_params()
            
            rec_rows.append({'id': rec_id, 'sync_rec_id': srec_id, 'cell_id': cell_id,
                'device_key': rec.device_id, 'stim_name': rec.meta['stim_name'],
                'clamp_mode': rec.clamp_mode, 'test_pulse_id': tp_id,
                'sample_rate': rec['primary'].sample_rate,
                'induction_freq': ind_freq, 'recovery_delay': recovery_delay})
            
            pulses = psa.pulses()
            if pulses[0][2] < 0:
                # test pulse
                tp = pulses[0]
            else:
                tp = (None, None)
            tp_rows.append({'id': tp_id, 'recording_id': rec_id,
                'pulse_start': tp[0], 'pulse_stop': tp[1],
                })
            
            if pre is None or rec.device_id == pre:
                pulse_n_id_map = {}
                for i,pulse in enumerate(pulses):
                    pulse_id = len(pulse_rows)
                    pulse_n_id_map[i] = pulse_id
                    pulse_key_n_id_map[(rec.device_id, i)] = pulse_id
                    pulse_rows.append({'id': pulse_id, 'recording_id': rec_id,
                        'pulse_number': i, 'onset_index': pulse[0],
                        'length': pulse[1]-pulse[0], 'amplitude': pulse[2],
                        'n_spikes': 0})
            
                spikes = psa.evoked_spikes()
                for sp in spikes:
                    sp_id = len(spike_rows)
                    pulse_id = pulse_n_id_map[sp['pulse_n']]
                    srow = {'id': sp_id, 'recording_id': rec_id, 'pulse_id': pulse_id}
                    pulse_rows[pulse_id]['n_spikes'] += 1
                    if sp['spike'] is not None:
                        srow.update(sp['spike'])
                    spike_rows.append(srow)
            
        mpa = MultiPatchSyncRecAnalyzer(srec)
        for pre_dev in srec.devices:
            if pre is not None and pre_dev != pre:
                continue
            
            for post_dev in srec.devices:
                if post is not None and post_dev != post:
                    continue
                
                responses = mpa.get_spike_responses(srec[pre_dev], srec[post_dev], align_to='pulse', require_spike=False)
                for resp in responses:
                    resp_id = len(response_rows)
                    bl_id = len(baseline_rows)
                    baseline_rows.append({'id': bl_id,
                        'recording_id': rec_key_id_map[post_dev],
                        'start_index': resp['baseline_start'],
                        'stop_index': resp['baseline_stop'],
                        'data': resp['baseline'].downsample(f=50000).data,
                        'value': float_mode(resp['baseline'].data),
                    })
                    response_rows.append({'id': resp_id, 
                        'recording_id': rec_key_id_map[post_dev],
                        'pulse_id': pulse_key_n_id_map[(pre_dev, resp['pulse_n'])],
                        'start_index': resp['rec_start'], 'stop_index': resp['rec_stop'],
                        'baseline_id': bl_id,
                        'data': resp['response'].downsample(f=50000).data,
                    })
    
    
def submit_site_mosaic(data):
    """Submit a site mosaic and information about cell labeling
    
        data = {
            'acquisition_site_uid': <unique site ID chosen by acquisition system>,
            'mosaic_file': <path to mosaic file>,
            'cells': [
                {'cell_id':, 'fill_fluorophore':, 'cre_labels': [], 'position': (x,y,z)},
            ],
        }
    
    """
    
    
def submit_biocytin_data(data):
    """Submit metadata related to biocytin image.
    
        data = {
            'cells': [
                {'cell_id':, 'biocytin_filled': bool},
            ],
        }
    
    """
    

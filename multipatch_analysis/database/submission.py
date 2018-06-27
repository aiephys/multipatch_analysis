from acq4.util.DataManager import getDirHandle
import os, re, json, yaml, shutil
from collections import OrderedDict
from datetime import datetime, timedelta
import numpy as np
import pyqtgraph as pg
from neuroanalysis.baseline import float_mode
from neuroanalysis.data import PatchClampRecording
from . import database as db
from .. import lims
from ..data import MultiPatchExperiment, MultiPatchProbe
from ..connection_detection import PulseStimAnalyzer, MultiPatchSyncRecAnalyzer, BaselineDistributor
from .. import config
from .. import constants
from .. import qc


class SliceSubmission(object):
    """Used to submit a new slice entry to the synphys DB.
    """
    message = "Generating slice DB entries"
    
    def __init__(self, slice_dir):
        self.slice_dir = slice_dir
        self.dh = getDirHandle(self.slice_dir)
        self._fields = None

    @property
    def fields(self):
        if self._fields is None:
            info = self.dh.info()
            
            # pull some metadata from LIMS
            sid = info['specimen_ID'].strip()
            limsdata = lims.specimen_info(sid)

            quality = info.get('slice quality', None)
            try:
                quality = int(quality)
            except Exception:
                quality = None

            # Interpret slice time
            slice_time = info.get('slice time', None)
            if slice_time is not None:
                m = re.match(r'((20\d\d)-(\d{1,2})-(\d{1,2}) )?(\d+):(\d+)', slice_time.strip())
                if m is not None:
                    _, year, mon, day, hh, mm = m.groups()
                    if year is None:
                        date = datetime.fromtimestamp(self.dh.parent().info('__timestamp__'))
                        slice_time = datetime(date.year, date.month, date.day, hh, mm)
                    else:
                        slice_time = datetime(year, mon, day, hh, mm)

            self._fields = {
                'acq_timestamp': datetime.fromtimestamp(info['__timestamp__']),
                'species': limsdata['organism'],
                'age': limsdata['age'],
                'sex': limsdata['sex'],
                'genotype': limsdata['genotype'],
                'orientation': limsdata['plane_of_section'],
                'surface': limsdata['exposed_surface'],
                'hemisphere': limsdata['hemisphere'],
                'quality': quality,
                'slice_time': slice_time,
                'slice_conditions': {},
                'lims_specimen_name': sid,
                'storage_path': self.dh.name(relativeTo=self.dh.parent().parent()),
                'submission_data': None,
            }
        return self._fields

    def check(self):
        warnings = []
        errors = []
        fields = self.fields
        
        ts = datetime.fromtimestamp(self.dh.info()['__timestamp__'])
        try:
            slice_entry = db.slice_from_timestamp(ts)
            update_fields = []
            for k,v in fields.items():
                v1 = getattr(slice_entry, k)
                if v1 != v:
                    update_fields.append((k, v1, v))
            if len(update_fields) > 0:
                update_msg = ['%s: %s=>%s'%x for x in update_fields]
                warnings.append("will overwrite slice metadata: %s" % (', '.join(update_msg)))
        except KeyError:
            pass
       
        if fields['surface'] not in ['right', 'left']:
            warnings.append("Warning: slice surface '%s' should have been 'right' or 'left'" % fields['surface'])
        
        return errors, warnings
        
    def submitted(self):
        slice_dir = self.dh
        ts = datetime.fromtimestamp(slice_dir.info()['__timestamp__'])
        try:
            slice_entry = db.slice_from_timestamp(ts)
            return True
        except KeyError:
            return False
        
    def summary(self):
        return {'slice': self.fields}
        
    def create(self):
        if len(self.check()[0]) > 0:
            raise Exception("Submission has errors; see SliceSubmission.check()")
        data = self.fields
        sl = db.Slice(**data)
        return sl
        
    def submit(self):
        session = db.Session()
        try:
            sl = self.create()
            session.add(sl)
            session.commit()
        finally:
            session.close()
        

class ExperimentDBSubmission(object):
    """Used to submit a new experiment entry to the synphys DB.
    
    This causes several tables to be populated: experiment, sync_rec, recording,
    stim_pulse, stim_spike, pulse_response, baseline
    """
    message = "Generating database entries"

    def __init__(self, expt):
        self.expt = expt
        self._fields = None

    def submitted(self):
        ts = self.expt.datetime
        try:
            expt_entry = db.experiment_from_timestamp(ts)
            return True
        except KeyError:
            return False

    def check(self):
        warnings = []
        errors = []
        expt = self.expt
        # info = self.dh.info()
        
        # slice_dir = self.dh.parent()
        
        # expt_dir = slice_dir.parent()
        # expt_info = expt_dir.info()
        # self._expt_info = expt_info

        try:
            temp = expt.target_temperature
        except:
            warnings.append("Experiment temperature '%s' is invalid." % self.expt.expt_info['temperature'])

        nwb_file = expt.nwb_file
        if not os.path.isfile(nwb_file):
            errors.append('Could not find NWB file "%s"' % nwb_file)
        
        expt_info = expt.expt_info

        self.fields = {
            'original_path': expt.original_path,
            'storage_path': expt.server_path,
            'ephys_file': os.path.relpath(expt.nwb_file, expt.path),
            'rig_name': expt.rig_name,
            'acq_timestamp': expt.datetime,
            'target_region': expt_info.get('region'),
            'internal': expt_info.get('internal'),
            'acsf': expt_info.get('solution'),
            'target_temperature': temp,
        }
        
        if self.submitted():
            errors.append("Experiment is already submitted.")

        return errors, warnings
    
    def summary(self):
        return {'database': 'might add some records..'}
        
    def create(self, session):
        err,warn = self.check()
        if len(err) > 0:
            raise Exception("Submission has errors:\n%s" % '\n'.join(err))

        # look up slice record in DB
        ts = self.expt.slice_timestamp
        slice_entry = db.slice_from_timestamp(ts, session=session)
        
        # Create entry in experiment table
        data = self.fields
        expt_entry = db.Experiment(**data)
        expt_entry.slice = slice_entry
        self.expt_entry = expt_entry
        session.add(expt_entry)

        # create pipette and cell entries
        elecs_by_ad_channel = {}
        cell_entries = {}
        for e_id, elec in self.expt.electrodes.items():
            elec_entry = db.Electrode(experiment=expt_entry, ext_id=elec.electrode_id, device_id=elec.device_id)
            for k in ['patch_status', 'start_time', 'stop_time',  
                      'initial_resistance', 'initial_current', 'pipette_offset',
                      'final_resistance', 'final_current']:
                if hasattr(elec, k):
                    setattr(elec_entry, k, getattr(elec, k))
            session.add(elec_entry)

            # store so recordings can reference this later on..
            elecs_by_ad_channel[elec.device_id] = elec_entry

            if elec.cell is not None:
                cell = elec.cell
                cell_entry = db.Cell(
                    electrode=elec_entry,
                    ext_id=cell.cell_id,
                    cre_type=cell.cre_type,
                    target_layer=cell.target_layer,
                    is_excitatory=cell.is_excitatory,
                    depth=cell.depth,
                    position=cell.position,
                )
                session.add(cell_entry)
                cell_entries[cell] = cell_entry

        # create pairs
        pairs_by_device_id = {}
        for i, pre_cell in self.expt.cells.items():
            for j, post_cell in self.expt.cells.items():
                if i == j:
                    continue
                # check to see if i,j is in manual connection calls
                # (do not use expt.connections, which excludes some connections based on QC)
                syn_calls = self.expt.connection_calls
                synapse = None if syn_calls is None else ((i, j) in syn_calls)
                gap_calls = self.expt.gap_calls
                electrical = None if gap_calls is None else ((i, j) in gap_calls)

                pre_cell_entry = cell_entries[pre_cell]
                post_cell_entry = cell_entries[post_cell]
                p1, p2 = pre_cell.position, post_cell.position
                if None in [p1, p2]:
                    distance = None
                else:
                    distance = np.linalg.norm(np.array(p1) - np.array(p2))
                
                pair_entry = db.Pair(
                    experiment=expt_entry,
                    pre_cell=pre_cell_entry,
                    post_cell=post_cell_entry,
                    synapse=synapse,
                    electrical=electrical,
                    n_ex_test_spikes=0,  # will be counted later
                    n_in_test_spikes=0,
                    distance=distance,
                )
                session.add(pair_entry)

                pre_id = pre_cell_entry.electrode.device_id
                post_id = post_cell_entry.electrode.device_id
                pairs_by_device_id[(pre_id, post_id)] = pair_entry

        # Load NWB file and create data entries
        self._load_nwb(session, expt_entry, elecs_by_ad_channel, pairs_by_device_id)

        return expt_entry

    def _load_nwb(self, session, expt_entry, elecs_by_ad_channel, pairs_by_device_id):
        nwb = self.expt.data
        
        for srec in nwb.contents:
            temp = srec.meta.get('temperature', None)
            srec_entry = db.SyncRec(ext_id=srec.key, experiment=expt_entry, temperature=temp)
            session.add(srec_entry)
            
            srec_has_mp_probes = False
            
            rec_entries = {}
            all_pulse_entries = {}
            for rec in srec.recordings:
                
                # import all recordings
                rec_entry = db.Recording(
                    sync_rec=srec_entry,
                    electrode=elecs_by_ad_channel[rec.device_id],  # should probably just skip if this causes KeyError?
                    start_time=rec.start_time,
                )
                session.add(rec_entry)
                rec_entries[rec.device_id] = rec_entry
                
                # import patch clamp recording information
                if not isinstance(rec, PatchClampRecording):
                    continue
                qc_pass = qc.recording_qc_pass(rec)
                pcrec_entry = db.PatchClampRecording(
                    recording=rec_entry,
                    clamp_mode=rec.clamp_mode,
                    patch_mode=rec.patch_mode,
                    stim_name=rec.stimulus.description,
                    baseline_potential=rec.baseline_potential,
                    baseline_current=rec.baseline_current,
                    baseline_rms_noise=rec.baseline_rms_noise,
                    qc_pass=qc_pass,
                )
                session.add(pcrec_entry)

                # import test pulse information
                tp = rec.nearest_test_pulse
                if tp is not None:
                    tp_entry = db.TestPulse(
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
                    pcrec_entry.nearest_test_pulse = tp_entry
                    
                # import information about STP protocol
                if not isinstance(rec, MultiPatchProbe):
                    continue
                srec_has_mp_probes = True
                psa = PulseStimAnalyzer.get(rec)
                ind_freq, rec_delay = psa.stim_params()
                mprec_entry = db.MultiPatchProbe(
                    patch_clamp_recording=pcrec_entry,
                    induction_frequency=ind_freq,
                    recovery_delay=rec_delay,
                )
                session.add(mprec_entry)
            
                # import presynaptic stim pulses
                pulses = psa.pulses()
                
                pulse_entries = {}
                all_pulse_entries[rec.device_id] = pulse_entries
                
                rec_tvals = rec['primary'].time_values

                for i,pulse in enumerate(pulses):
                    # Record information about all pulses, including test pulse.
                    t0 = rec_tvals[pulse[0]]
                    t1 = rec_tvals[pulse[1]]
                    data_start = max(0, t0 - 10e-3)
                    data_stop = t0 + 10e-3
                    pulse_entry = db.StimPulse(
                        recording=rec_entry,
                        pulse_number=i,
                        onset_time=t0,
                        amplitude=pulse[2],
                        duration=t1-t0,
                        data=rec['primary'].time_slice(data_start, data_stop).resample(sample_rate=20000).data,
                        data_start_time=data_start,
                    )
                    session.add(pulse_entry)
                    pulse_entries[i] = pulse_entry
                    

                # import presynaptic evoked spikes
                # For now, we only detect up to 1 spike per pulse, but eventually
                # this may be adapted for more.
                spikes = psa.evoked_spikes()
                for i,sp in enumerate(spikes):
                    pulse = pulse_entries[sp['pulse_n']]
                    if sp['spike'] is not None:
                        spinfo = sp['spike']
                        extra = {
                            'peak_time': rec_tvals[spinfo['peak_index']],
                            'max_dvdt_time': rec_tvals[spinfo['rise_index']],
                            'max_dvdt': spinfo['max_dvdt'],
                        }
                        if 'peak_diff' in spinfo:
                            extra['peak_diff'] = spinfo['peak_diff']
                        if 'peak_value' in spinfo:
                            extra['peak_value'] = spinfo['peak_value']
                        
                        pulse.n_spikes = 1
                    else:
                        extra = {}
                        pulse.n_spikes = 0
                    
                    spike_entry = db.StimSpike(
                        pulse=pulse,
                        **extra
                    )
                    session.add(spike_entry)
                    pulse.first_spike = spike_entry
            
            if not srec_has_mp_probes:
                continue
            
            # import postsynaptic responses
            mpa = MultiPatchSyncRecAnalyzer(srec)
            for pre_dev in srec.devices:
                for post_dev in srec.devices:
                    if pre_dev == post_dev:
                        continue

                    # get all responses, regardless of the presence of a spike
                    responses = mpa.get_spike_responses(srec[pre_dev], srec[post_dev], align_to='pulse', require_spike=False)
                    post_tvals = srec[post_dev]['primary'].time_values
                    for resp in responses:
                        # base_entry = db.Baseline(
                        #     recording=rec_entries[post_dev],
                        #     start_index=resp['baseline_start'],
                        #     stop_index=resp['baseline_stop'],
                        #     data=resp['baseline'].resample(sample_rate=20000).data,
                        #     mode=float_mode(resp['baseline'].data),
                        # )
                        # session.add(base_entry)
                        pair_entry = pairs_by_device_id.get((pre_dev, post_dev), None)
                        if pair_entry is None:
                            continue  # no data for one or both channels
                        if resp['ex_qc_pass']:
                            pair_entry.n_ex_test_spikes += 1
                        if resp['in_qc_pass']:
                            pair_entry.n_in_test_spikes += 1
                        resp_entry = db.PulseResponse(
                            recording=rec_entries[post_dev],
                            stim_pulse=all_pulse_entries[pre_dev][resp['pulse_n']],
                            pair=pair_entry,
                            # baseline=base_entry,
                            start_time=post_tvals[resp['rec_start']],
                            data=resp['response'].resample(sample_rate=20000).data,
                            ex_qc_pass=resp['ex_qc_pass'],
                            in_qc_pass=resp['in_qc_pass'],
                        )
                        session.add(resp_entry)
                        
            # generate up to 20 baseline snippets for each recording
            for dev in srec.devices:
                rec = srec[dev]
                rec_tvals = rec['primary'].time_values
                dist = BaselineDistributor.get(rec)
                for i in range(20):
                    base = dist.get_baseline_chunk(20e-3)
                    if base is None:
                        # all out!
                        break
                    start, stop = base
                    data = rec['primary'][start:stop].resample(sample_rate=20000).data

                    ex_qc_pass, in_qc_pass = qc.pulse_response_qc_pass(rec, [start, stop], None, [])

                    base_entry = db.Baseline(
                        recording=rec_entries[dev],
                        start_time=rec_tvals[start],
                        data=data,
                        mode=float_mode(data),
                        ex_qc_pass=ex_qc_pass,
                        in_qc_pass=in_qc_pass,
                    )
                    session.add(base_entry)
            
        
    def submit(self):
        session = db.Session()
        try:
            exp = self.create(session)
            session.commit()
        except:
            session.rollback()
            raise
        finally:
            session.close()

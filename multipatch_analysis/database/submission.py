from acq4.util.DataManager import getDirHandle
import os, re, json, yaml, shutil
from collections import OrderedDict
from datetime import datetime, timedelta
import pyqtgraph as pg
from neuroanalysis.baseline import float_mode
from neuroanalysis.data import PatchClampRecording
from . import database as db
from .. import lims
from ..data import MultiPatchExperiment, MultiPatchProbe
from ..connection_detection import PulseStimAnalyzer, MultiPatchSyncRecAnalyzer, BaselineDistributor
from .. import config
from .. import constants


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
                m = re.match(r'((20\d\d)-(\d{1,2})-(\d{1,2})\s+)?(\d+):(\d+)', slice_time.strip())
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
                )
                session.add(cell_entry)
                cell_entries[cell] = cell_entry

        # create pairs
        pairs_by_device_id = {}
        for i, pre_cell in self.expt.cells.items():
            for j, post_cell in self.expt.cells.items():
                if i == j:
                    continue
                synapse = (i, j) in self.expt.connections
                pre_cell_entry = cell_entries[pre_cell]
                post_cell_entry = cell_entries[post_cell]
                pair_entry = db.Pair(
                    experiment=expt_entry,
                    pre_cell=pre_cell_entry,
                    post_cell=post_cell_entry,
                    synapse=synapse,
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
            srec_entry = db.SyncRec(sync_rec_key=srec.key, experiment=expt_entry, temperature=temp)
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
                        resp_entry = db.PulseResponse(
                            recording=rec_entries[post_dev],
                            stim_pulse=all_pulse_entries[pre_dev][resp['pulse_n']],
                            pair=pairs_by_device_id[(pre_dev, post_dev)],
                            # baseline=base_entry,
                            start_time=post_tvals[resp['rec_start']],
                            data=resp['response'].resample(sample_rate=20000).data,
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
                    start, stop = base.source_indices
                    base_entry = db.Baseline(
                        recording=rec_entries[dev],
                        start_time=rec_tvals[start],
                        data=base.resample(sample_rate=20000).data,
                        mode=float_mode(base.data),
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


class PipetteDBSubmission(object):
    """Generates electrode / cell / pair entries in the synphys DB
    """
    message = "Generating pipette entries in DB"
    
    def __init__(self, site_dh):
        self.site_dh = site_dh
        
    def check(self):
        errors = []
        info = self.site_dh.info()
        if 'pipettes' not in info:
            errors.append("No pipette information found in site metadata.")
        return errors, []
        
    def summary(self):
        return {'pipettes': None}
        
    def submit(self):
        ts = datetime.fromtimestamp(self.site_dir.info()['__timestamp__'])
        expt_entry = db.experiment_from_timestamp(ts)


class SiteMosaicSubmission(object):
    """Submit site mosaic data to multiple locations:
    
    * Copy site mosaic file to server
    * Add cell position/label information to DB
    * Update LIMS json
    """
    message = "Submitting site mosaic"
    
    def __init__(self, site_dh):
        self.site_dh = site_dh
        
    def check(self):
        return [], []
        
    def summary(self):
        return {'mosaic': None}
        
    def submit(self):
        pass


class ExperimentSubmission(object):
    def __init__(self, site_dh, files, pipettes):
        """Handles the entire process of submitting an experiment, including:
        
        1. Storing submission metadata into ACQ4 index files
        2. Copying all files to synphys data storage
        3. Submitting combined NWB + metadata to LIMS
        
        Parameters
        ----------
        dh : DirHandle
            Handle to the ACQ4 "site" folder to be submitted
        files : list
            Describes all files to be considered part of this experiment.
            Each item in the list is a dict containing:
               
                * path: file path relative to slice folder
                * category: the purpose of this file
        pipettes : list
            Describes all electrodes used during this experiment. Each item in
            the list is a dict containing:
            
                * id: the ID of this pipette (usually 1-8)
                * status: 'GOhm seal', 'Low seal', 'No seal', ...
                * got_cell: bool, whether a cell was recorded from
                * channel: AD channel number for this pipette
                * start: starting datetime
                * stop: ending datetime
        """
        self.site_dh = site_dh
        self.files = files
        self.pipettes = pipettes
    
        nwb_file = [f['path'] for f in files if f['category'] == 'MIES physiology'][0]
        self.nwb_file = site_dh.parent()[nwb_file]
    
        self.stages = [
            ExperimentMetadataSubmission(site_dh, files, pipettes),
            #RawDataSubmission(site_dh),
            #LIMSSubmission(site_dh, files),
            #LIMSSubmission(site_dh, files, _override_spec_name="Ntsr1-Cre_GN220;Ai14-349905.03.06"),
        ]
        
    def check(self):
        errors = []
        warnings = []

        # sort files by category for counting
        categories = {}
        for f in self.files:
            categories.setdefault(f['category'], []).append(f['path'])
        n_files = {k:len(v) for k,v in categories.items()}
        
        # Make sure exactly one NWB file was selected
        n_mp = n_files.get('MIES physiology', 0)
        if n_mp != 1:
            errors.append('selected %d NWB files (must have exactly 1)' % n_mp)
        
        # Make sure exactly one MP log file was selected
        n_mpl = n_files.get('Multipatch log', 0)
        if n_mpl != 1:
            errors.append('selected %d multipatch log files (must have exactly 1)' % n_files['Multipatch log'])

        # Check that we have at least 1 slice anatomy image, and warn if we have fewer
        # than 4
        n_sa = n_files.get('slice anatomy', 0)
        if n_sa < 1:
            errors.append('no slice anatomy images found.')
        elif n_sa < 4:
            warnings.append('only %d slice anatomy images found.' % n_sa)

        # warn if there are no slice quality images
        n_sq = n_files.get('slice quality stack', 0)
        if n_sq < 1:
            warnings.append('no slice quality stack found.')
            
        # todo: check for images covering all fluorescent reporters

        # todo: make sure all images are associated with a light source
        
        # allow all submission stages to make their own checks
        for stage in self.stages:
            e,w = stage.check()
            errors.extend(e)
            warnings.extend(w)
            
        return errors, warnings

    def summary(self):
        summ = OrderedDict()
        for stage in self.stages:
            summ.update(stage.summary())
        return summ
    
    def submit(self):
        for stage in self.stages:
            assert len(stage.check()[0]) == 0
        with pg.ProgressDialog("Submitting experiment..", maximum=len(self.stages), nested=True) as dlg:
            for stage in self.stages:
                stage.submit()
                dlg += 1
                if dlg.wasCanceled():
                    break


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
    

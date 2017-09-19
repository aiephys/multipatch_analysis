from collections import OrderedDict
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
            warnings.append("Warning: slice surface '%s' should have been 'right' or 'left'" % surface)
        
        return errors, warnings
        
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
    
    
class ExperimentDBSubmission(object):
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
            
            expt_dir = slice_dir.parent()
            expt_info = expt_dir.info()

            temp = expt_info.get('temperature')
            if temp is not None:
                temp = float(temp.rstrip(' C'))

            self._fields = {
                'original_path': '%s:%s' % (config.rig_name, self.dh.name()),
                'acq_timestamp': datetime.fromtimestamp(info['__timestamp__']),
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
    
    def summary(self):
        return {'database': 'might add some records..'}
        
    def create(self, session):
        if len(self.check()[0]) > 0:
            raise Exception("Submission has errors; see SiteSubmission.check()")

        # look up slice record in DB
        slice_dir = self.dh.parent()
        ts = datetime.fromtimestamp(slice_dir.info()['__timestamp__'])
        slice_entry = db.slice_from_timestamp(ts)

        
        # Create entry in experiment table
        data = self.fields
        expt = db.Experiment(**data)
        expt.slice = slice_entry
        self.expt_entry = expt
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


class ExperimentMetadataSubmission(object):
    """Handles storing submission metadata to ACQ4 index files.
    
    (see ExperimentSubmission for details)
    """
    def __init__(self, site_dh, files, electrodes):
        self.site_dh = site_dh
        self.files = files
        self.electrodes = electrodes
        
    def check(self):
        """
        * error if any files do not exist
        * warning if existing files are not mentioned in file list
        * warning if any metadata would be overwritten with a new value
        """
        errors = []
        warnings = []
        
        slice_dh = self.site_dh.parent()
        for file_info in self.files:
            fh = slice_dh[file_info['path']]
            if not fh.isFile():
                errors.append("file %s does not exist" % fh.name())
                continue
            info = fh.info()
            if 'category' in info and info['category'] != file_info['category']:
                warnings.append("category for file %s will change %s => %s" % 
                                (info['category'], file_info['category']))
        
        all_files = [f['path'] for f in self.files]
        for dh in [slice_dh, self.site_dh]:
            for fname in dh.ls():
                fh = dh[fname]
                relname = fh.name(relativeTo=slice_dh)
                if fh.isFile() and relname not in all_files:
                    warnings.append("file %s is not mentioned in metadata." % relname)
        
        site_info = self.site_dh.info()
        if 'electrodes' in site_info and site_info['electrodes'] != self.electrodes:
            warnings.append('site electrode metadata will be overwritten: %s => %s' % 
                            (site_info['electrodes'], self.electrodes))
        
        return errors, warnings

    def summary(self):
        summ = OrderedDict()
        summ['file categories'] = self.files
        summ['electrodes'] = self.electrodes
        
        return {'metadata': summ}
    
    def submit(self):
        slice_dh = self.site_dh.parent()
        for file_info in self.files:
            fh = slice_dh[file_info['path']]
            fh.setInfo(category=file_info['category'])
        
        self.site_dh.setInfo(electrodes=self.electrodes)
        
            
        
        
        


class ExperimentSubmission(object):
    def __init__(self, site_dh, files, electrodes):
        """Handles the entire process of submitting an experiment, including:
        
        1. Storing submission metadata into ACQ4 index files
        2. Submitting slice metadata to synphys DB
        3. Submitting initial experiment to synphys DB
        4. Submitting electrode / cell info to synphys DB
        5. Copying all files to synphys data storage
        6. Generating NWB file
        7. Submitting NWB + metadata to LIMS
        
        Parameters
        ----------
        dh : DirHandle
            Handle to the ACQ4 "site" folder to be submitted
        files : list
            Describes all files to be considered part of this experiment.
            Each item in the list is a dict containing:
               
                * path: file path relative to slice folder
                * category: the purpose of this file
        electrodes : list
            Describes all electrodes used during this experiment. Each item in
            the list is a dict containing:
            
                * id: the ID of this electrode (usually 1-8)
                * status: 'GOhm seal', 'Low seal', 'No seal', ...
                * got_cell: bool, whether a cell was recorded from
                * channel: AD channel number for this electrode
                * start: starting datetime
                * stop: ending datetime
        """
        self.site_dh = site_dh
        self.files = files
        self.electrodes = electrodes
    
        nwb_file = [f['path'] for f in files if f['category'] == 'MIES physiology'][0]
        self.nwb_file = site_dh.parent()[nwb_file]
    
        self.stages = [
            ExperimentMetadataSubmission(site_dh, files, electrodes),
            SliceSubmission(site_dh.parent()),
            ExperimentDBSubmission(site_dh, self.nwb_file),
        ]
        
    def check(self):
        errors = []
        warnings = []
        
        # Make sure only one NWB file was selected
        nwb_files = [f['path'] for f in self.files if f['category'] == 'MIES physiology']
        if len(nwb_files) != 1:
            errors.append('selected %d NWB files (must have exactly 1)' % len(nwb_files))
        
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
        for stage in self.stages:
            stage.submit()


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
    

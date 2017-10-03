import os, re, json, shutil
from collections import OrderedDict
from datetime import datetime
import database as db
from neuroanalysis.baseline import float_mode
from neuroanalysis.data import PatchClampRecording
import lims
from data import MultiPatchExperiment, MultiPatchProbe
from connection_detection import PulseStimAnalyzer, MultiPatchSyncRecAnalyzer
import config


class SliceSubmission(object):
    """Used to submit a new slice entry to the synphys DB.
    """
    message = "Generating slice DB entries"
    
    def __init__(self, dh):
        self.dh = dh
        
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
                'genotype': limsdata['genotype'],
                'orientation': limsdata['plane_of_section'],
                'surface': limsdata['exposed_surface'],
                'hemisphere': limsdata['hemisphere'],
                'quality': quality,
                'slice_time': slice_time,
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
    message = "Generating database entries"

    def __init__(self, dh, nwb_file):
        self.dh = dh
        self.nwb_file = nwb_file
        self._fields = None

    def submitted(self):
        site_dir = self.dh
        ts = datetime.fromtimestamp(site_dir.info()['__timestamp__'])
        try:
            expt_entry = db.experiment_from_timestamp(ts)
            return True
        except KeyError:
            return False

    def check(self):
        warnings = []
        errors = []
        
        info = self.dh.info()
        
        slice_dir = self.dh.parent()
        
        expt_dir = slice_dir.parent()
        expt_info = expt_dir.info()
        self._expt_info = expt_info

        temp = expt_info.get('temperature')
        if temp is not None:
            temp = temp.lower().rstrip(' c').strip()
            if temp == 'rt':
                temp = 22.0
            else:
                try:
                    temp = float(temp)
                except Exception:
                    temp = None
                    errors.append("Experiment temperature '%s' is invalid." % self._expt_info['temperature'])
        
        self.fields = {
            'original_path': '%s:%s' % (config.rig_name, self.dh.name()),
            'acq_timestamp': datetime.fromtimestamp(info['__timestamp__']),
            'target_region': expt_info.get('region'),
            'internal': expt_info.get('internal'),
            'acsf': expt_info.get('acsf'),
            'target_temperature': temp,
            
        }
                
        return errors, warnings
    
    def summary(self):
        return {'database': 'might add some records..'}
        
    def create(self, session):
        err,warn = self.check()
        if len(err) > 0:
            raise Exception("Submission has errors:\n%s" % '\n'.join(err))

        # look up slice record in DB
        slice_dir = self.dh.parent()
        ts = datetime.fromtimestamp(slice_dir.info()['__timestamp__'])
        slice_entry = db.slice_from_timestamp(ts, session=session)
        
        # Create entry in experiment table
        data = self.fields
        expt = db.Experiment(**data)
        expt.slice = slice_entry
        self.expt_entry = expt
        session.add(expt)
        
        # Load NWB file and create data entries
        nwb = MultiPatchExperiment(self.nwb_file.name())
        
        for srec in nwb.contents:
            temp = srec.meta.get('temperature', None)
            srec_entry = db.SyncRec(sync_rec_key=srec.key, experiment=expt, temperature=temp)
            session.add(srec_entry)
            
            srec_has_mp_probes = False
            
            rec_entries = {}
            all_pulse_entries = {}
            for rec in srec.recordings:
                
                # import all recordings
                rec_entry = db.Recording(
                    sync_rec=srec_entry,
                    device_key=rec.device_id, 
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
                
                for i,pulse in enumerate(pulses):
                    # Record information about all pulses, including test pulse.
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
            
            if not srec_has_mp_probes:
                continue
            
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
                            data=resp['baseline'].resample(sample_rate=20000).data,
                            mode=float_mode(resp['baseline'].data),
                        )
                        session.add(base_entry)
                        resp_entry = db.PulseResponse(
                            recording=rec_entries[post_dev],
                            stim_pulse=all_pulse_entries[pre_dev][resp['pulse_n']],
                            baseline=base_entry,
                            start_index=resp['rec_start'],
                            stop_index=resp['rec_stop'],
                            data=resp['response'].resample(sample_rate=20000).data,
                        )
                        session.add(resp_entry)
                        
        return expt
        
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


class ExperimentMetadataSubmission(object):
    """Handles storing submission metadata back to ACQ4 index files.
    
    This makes it possible to automatically re-submit without user intervention.
    
    (see ExperimentSubmission for details)
    """
    message = "Copying submission metadata back to index files"
    
    def __init__(self, site_dh, files, pipettes):
        self.site_dh = site_dh
        self.files = files
        self.pipettes = pipettes
        
    def check(self):
        """
        * error if any files do not exist
        * warning if existing files are not mentioned in file list
        * warning if any metadata would be overwritten with a new value
        """
        errors = []
        warnings = []
        
        # Do all categorized files actually exist?
        # Was a file categorized differently on a previous submission?
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
        
        # Did we forget to categorize a file?
        all_files = [f['path'] for f in self.files]
        for dh in [slice_dh, self.site_dh]:
            for fname in dh.ls():
                fh = dh[fname]
                relname = fh.name(relativeTo=slice_dh)
                if fh.isFile() and relname not in all_files:
                    warnings.append("file %s is not mentioned in metadata." % relname)
        
        # Do we already have pipette information submitted?
        site_info = self.site_dh.info()
        if 'pipettes' in site_info and site_info['pipettes'] != self.pipettes:
            warnings.append('site pipette metadata will be overwritten: %s => %s' % 
                            (site_info['pipettes'], self.pipettes))
        
        return errors, warnings

    def summary(self):
        summ = OrderedDict()
        summ['file categories'] = self.files
        summ['pipettes'] = self.pipettes
        
        return {'metadata': summ}
    
    def submit(self):
        slice_dh = self.site_dh.parent()
        expt_dh = slice_dh.parent()
        expt_dh.setInfo(rig_name=config.rig_name)
        for file_info in self.files:
            if file_info['category'] == 'ignore':
                continue
            fh = slice_dh[file_info['path']]
            fh.setInfo(category=file_info['category'])
        
        self.site_dh.setInfo(pipettes=self.pipettes)


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


class RawDataSubmission(object):
    """Copies all raw data to a central server.
    """
    message = "Copying data to server"
    
    def __init__(self, site_dh):
        self.site_dh = site_dh
        
    def check(self):
        return [], []
        
    def summary(self):
        return {'raw_data': None}
        
    def submit(self):
        site_dh = self.site_dh
        slice_dh = site_dh.parent()
        expt_dh = slice_dh.parent()
        
        print("Copying raw data to server..")
        
        # Decide how the top-level directory will be named on the remote server
        # (it may already be there from a previous slice/site, or the current
        # name may already be taken by another rig.)
        server_path = config.synphys_data
        expt_name = expt_dh.info().get('synphys_server_path', None)
        if expt_name is None:
            expt_base_name = expt_dh.shortName().split('_')[0]
            expt_dirs = set(os.listdir(server_path))
            i = 0
            while True:
                expt_name = expt_base_name + '_%03d' % i
                if expt_name not in expt_dirs:
                    break
                i += 1
        server_expt_path = os.path.join(server_path, expt_name)
        
        # Leave an artifact
        expt_dh.setInfo(synphys_server_path=server_expt_path)

        self.server_path = server_expt_path
        print("   using server path: %s" % server_path)
        self._sync_paths(expt_dh.name(), server_expt_path)
        
        # Copy slice files if needed
        server_slice_path = os.path.join(server_expt_path, slice_dh.shortName())
        self._sync_paths(slice_dh.name(), server_slice_path)

        # Copy site files if needed
        server_site_path = os.path.join(server_slice_path, site_dh.shortName())
        self._sync_paths(site_dh.name(), server_site_path)

        # Leave an artifact indicating this copy finished
        site_dh.setInfo(synphys_server_path=server_site_path)
        
        print("Done.")
        
    def _sync_paths(self, source, target):
        """Non-recursive directory sync.
        """
        if not os.path.isdir(target):
            os.mkdir(target)
        for fname in os.listdir(source):
            src_path = os.path.join(source, fname)
            if os.path.isfile(src_path):
                dst_path = os.path.join(target, fname)
                if os.path.isfile(dst_path):
                    # Todo: add better checking here..
                    continue
                shutil.copyfile(src_path, dst_path)
                print("    copy %s => %s" % (src_path, dst_path))
        

class LIMSSubmission(object):
    """Packages all raw data into a single NWB/JSON pair and submits to LIMS.
    """
    message = "Submitting data to LIMS"
    
    def __init__(self, site_dh, files, _override_spec_name=None):
        self.site_dh = site_dh
        self.files = files
        
        # Only used for testing -- we can override the specimen name
        # in order to force the data to be uploaded to a specific specimen
        # that will not be used for real data.
        self.spec_name = _override_spec_name
        
    def check(self):
        errors = []
        warnings = []

        site_info = self.site_dh.info()
        slice_info = self.site_dh.parent().info()
        self.meta = {
            'acq_timestamp': site_info['__timestamp__'],
            'source_path': '%s:%s' % (config.rig_name, self.site_dh.name()),
        }
        
        if self.spec_name is None:
            self.spec_name = self.site_dh.parent().info()['specimen_ID'].strip()
        
        try:
            sid = lims.specimen_id_from_name(self.spec_name)
        except ValueError as err:
            errors.append(err.message)
            # bail out here; downstream checks will fail.
            return errors, warnings
        
        # LIMS upload will fail if the specimen has not been given an ephys roi plan.
        roi_plans = lims.specimen_ephys_roi_plans(self.spec_name)
        lims_edit_href = '<a href="http://lims2/specimens/{sid}/edit">http://lims2/specimens/{sid}/edit</a>'.format(sid=sid)
        if len(roi_plans) == 0:
            errors.append('Specimen has no ephys roi plan. Edit:' + lims_edit_href)
        elif len(roi_plans) == 1 and roi_plans[0]['name'] != 'Synaptic Physiology ROI Plan':
            errors.append('Specimen has wrong ephys roi plan '
                '(expected "Synaptic Physiology ROI Plan"). Edit:' + lims_edit_href)
        elif len(roi_plans) > 1:
            errors.append('Specimen has multiple ephys roi plans '
                '(expected 1). Edit:' + lims_edit_href)

        # Check LIMS specimen has flipped field set and a sensible-looking 
        # histology well name
        limsdata = lims.specimen_info(self.spec_name)
        if limsdata['flipped'] not in (True, False):
            errors.append("Must set flipped field for this specimen: %s" % lims_edit_href)
        # histology well name should look something like "multi_170911_21_A01"
        hist_well = limsdata['histology_well_name']
        if hist_well is None:
            errors.append("Missing histology well name in LIMS: %s" % lims_edit_href)        
        else:
            m = re.match(r'multi_(\d{2})(\d{2})(\d{2})_(2[1-9])_([A-C]0[1-4])', hist_well)
            if m is None:
                errors.append("Histology well name appears to be incorrectly formatted: %s" % lims_edit_href)        
            else:
                yy, mm, dd, plate_n, well = m.groups()[:3]
                plate_date = datetime(2000+int(yy), int(mm), int(dd))
                site_date = datetime.fromtimestamp(self.meta['acq_timestamp'])
                # find the most recent Monday
                expected_plate_date = site_date - timedelta(days=site_date.weekday())
                if abs((expected_plate_date - plate_date).total_seconds()) > 7*24*3600:
                    # error if more than a week out of sync
                    errors.append("Histology well date is %d%d%d; expected %s: %s" % (yy, mm, dd, expected_plate_date.strftime('%y%m%d'), lims_edit_href))
                if expected_plate_date.date() != plate_date:
                    # warning if the date is not exactly as expected
                    warnings.append("Histology well date is %d%d%d; expected %s: %s" % (yy, mm, dd, expected_plate_date.strftime('%y%m%d'), lims_edit_href))
                    
                if int(plate_n) > 24:
                    warnings.append("Histology plate number %s is probably too high. %s" % (plate_n, lims_edit_href))
                elif int(plate_n) > 22:
                    warnings.append("Histology plate number %s might be too high? %s" % (plate_n, lims_edit_href))

        # Check carousel ID matches the one in LIMS
        cw_name = limsdata['carousel_well_name']
        if cw_name is None or len(cw_name.strip()) == 0:
            errors.append('No LIMS carousel well name for this specimen!')
        else:
            if cw_name != slice_info['carousel_well_ID'].strip():
                errors.append('LIMS carousel well name "%s" does not match ACQ4 carousel_well_ID "%s"' 
                    % (cw_name, slice_info['carousel_well_ID']))
            
        # If histology well name was recorded in ACQ4, make sure it matches the one in LIMS
        acq_plate_well = slice_info.get('plate_well_ID', None)
        if acq_plate_well is not None and acq_plate_well.strip() != hist_well:
            errors.append('LIMS histology well name "%s" does not match ACQ4 plate_well_ID "%s"' 
                    % (hist_well, acq_plate_well))

        # make sure we have one nwb file to submit
        source_nwb_files = [f['path'] for f in self.files if f['category'] == 'MIES physiology']
        if len(source_nwb_files) == 0:
            errors.append("%d NWB files specified (should be 1)" % len(source_nwb_files))
        self.source_nwb = source_nwb_files[0]
        
        return errors, warnings
        
    def summary(self):
        return {'lims': None}
        
    def submit(self):
        import nwb_packaging
        
        assert len(self.check()[0]) == 0
        
        # Generate combined NWB file
        print("Generating combined NWB file..")
        slice_dh = self.site_dh.parent()
        file_paths = [slice_dh[f['path']].name() for f in self.files if f['path'] is not self.source_nwb]
        source_nwb = slice_dh[self.source_nwb].name()
        combined_nwb = nwb_packaging.buildCombinedNWB(source_nwb, file_paths)
        assert os.path.isfile(combined_nwb)
        
        # Generate json metadata file
        json_file = os.path.splitext(combined_nwb)[0] + '.json'
        json.dump(self.meta, open(json_file, 'wb'))
        print("Generated files:")
        print("    " + combined_nwb)
        print("    " + json_file)
        
        # submit to LIMS
        print("Submitting to LIMS..")
        self.lims_incoming_files = lims.submit_expt(self.spec_name, str(self.meta['acq_timestamp']), combined_nwb, json_file)
        
        # delete submission files
        os.remove(json_file)
        os.remove(combined_nwb)
        
        # Write an artifact file describing the upload data.
        # Ideally, we can use this to easily repeat the LIMS upload
        # in the future.
        artifact = {
            'files': self.files,
            'meta': self.meta,
            'specimen': self.spec_name,
            'submitted_files': self.lims_incoming_files,
        }
        artifact_file = self.site_dh['lims_upload_manifest.json'].name()
        json.dump(artifact, open(artifact_file, 'wb'))
        print("Complete; left artifact in %s" % artifact_file)


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
            RawDataSubmission(site_dh),
            LIMSSubmission(site_dh, files),
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
    

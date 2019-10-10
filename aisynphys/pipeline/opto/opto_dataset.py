from collections import OrderedDict
import os
from ... import config, qc #, synphys_cache, lims
from ...util import timestamp_to_datetime
from ..pipeline_module import DatabasePipelineModule
from .opto_experiment import OptoExperimentPipelineModule
from neuroanalysis.data.experiment import Experiment
from neuroanalysis.data.libraries import opto
#from neuroanalysis.stimuli import find_square_pulses
from ...data import PulseStimAnalyzer #Experiment, MultiPatchDataset, MultiPatchProbe, PulseStimAnalyzer, MultiPatchSyncRecAnalyzer, BaselineDistributor
from neuroanalysis.data import PatchClampRecording
from optoanalysis.optoadapter import OptoRecording
import optoanalysis.power_calibration as power_cal


class OptoDatasetPipelineModule(DatabasePipelineModule):

    name='opto_dataset'
    dependencies = [OptoExperimentPipelineModule]
    table_group = [
        'sync_rec', 
        'recording', 
        'patch_clamp_recording', 
        'test_pulse', 
        'stim_pulse', 
        'pulse_response',
        'baseline'
        ]

    # datasets are large and NWB access leaks memory -- this is probably true for optoanalysis too
    # when running parallel, each child process may run only one job before being killed
    maxtasksperchild = 1 

    @classmethod
    def create_db_entries(cls, job, session):
        db = job['database']
        job_id = job['job_id']

        # Load experiment from DB
        expt_entry = db.experiment_from_ext_id(job_id, session=session)
        elecs_by_ad_channel = {elec.device_id:elec for elec in expt_entry.electrodes}

        # load NWB file
        path = os.path.join(config.synphys_data, expt_entry.storage_path)
        expt = Experiment(site_path=path, loading_library=opto)
        nwb = expt.data
        stim_log = expt.library.load_stimulation_log(expt)

        #if expt.ephys_file is not None:
        #    raise Exception('stop')
        #else:
        #    return

        # Load all data from NWB into DB
        for srec in nwb.contents:
            temp = srec.meta.get('temperature', None)
            srec_entry = db.SyncRec(ext_id=srec.key, experiment=expt_entry, temperature=temp)
            session.add(srec_entry)

            rec_entries = {}
            all_pulse_entries = {}
            for rec in srec.recordings:
                
                # import all recordings
                electrode_entry = elecs_by_ad_channel.get(rec.device_id, None)
                #if electrode_entry is not None:
                rec_entry = db.Recording(
                    sync_rec=srec_entry,
                    electrode=electrode_entry,
                    start_time=rec.start_time,
                    device_name=str(rec.device_id)
                )
                session.add(rec_entry)
                rec_entries[rec.device_id] = rec_entry
                # else:
                #     opto_rec_entry = db.OpticalRecoding(
                #         sync_rec=srec_entry,
                #         device_name=str(rec.device_id),
                #         start_time=rec.start_time
                #     )
                #     session.add(opto_rec_entry)

                # import patch clamp recording information
                if isinstance(rec, PatchClampRecording):
                    qc_pass, qc_failures = qc.recording_qc_pass(rec)
                    pcrec_entry = db.PatchClampRecording(
                        recording=rec_entry,
                        clamp_mode=rec.clamp_mode,
                        patch_mode=rec.patch_mode,
                        stim_name=rec.stimulus.description,
                        baseline_potential=rec.baseline_potential,
                        baseline_current=rec.baseline_current,
                        baseline_rms_noise=rec.baseline_rms_noise,
                        qc_pass=qc_pass,
                        meta=None if len(qc_failures) == 0 else {'qc_failures': qc_failures},
                    )
                    session.add(pcrec_entry)

                    # import test pulse information
                    tp = rec.nearest_test_pulse
                    if tp is not None:
                        indices = tp.indices or [None, None]
                        tp_entry = db.TestPulse(
                            electrode=electrode_entry,
                            recording=rec_entry,
                            start_index=indices[0],
                            stop_index=indices[1],
                            baseline_current=tp.baseline_current,
                            baseline_potential=tp.baseline_potential,
                            access_resistance=tp.access_resistance,
                            input_resistance=tp.input_resistance,
                            capacitance=tp.capacitance,
                            time_constant=tp.time_constant,
                        )
                        session.add(tp_entry)
                        pcrec_entry.nearest_test_pulse = tp_entry

                    psa = PulseStimAnalyzer.get(rec)
                    pulses = psa.pulse_chunks()
                    pulse_entries = {}
                    all_pulse_entries[rec.device_id] = pulse_entries

                    for i,pulse in enumerate(pulses):
                        # Record information about all pulses, including test pulse.
                        t0, t1 = pulse.meta['pulse_edges']
                        resampled = pulse['primary'].resample(sample_rate=20000)
                        pulse_entry = db.StimPulse(
                            recording=rec_entry,
                            pulse_number=pulse.meta['pulse_n'],
                            onset_time=t0,
                            amplitude=pulse.meta['pulse_amplitude'],
                            duration=t1-t0,
                            data=resampled.data,
                            data_start_time=resampled.t0,
                            cell=electrode_entry.cell,
                            device_name=str(rec.device_id)
                        )
                        session.add(pulse_entry)
                        pulse_entries[pulse.meta['pulse_n']] = pulse_entry

                elif isinstance(rec, OptoRecording) and (rec.device_id=='Fidelity' or 'AD6' in rec.device_id): 
                    #psa = PulseStimAnalyzer.get(rec)
                    #psa.pulses(channel='reporter')
                    #pulses = psa.pulse_chunks
                    #pulses = find_square_pulses(rec['reporter'])
                    print('adding OptoRecording')
                    stim_num = rec.meta['notebook']['USER_stim_num']
                    stim = stim_log[str(int(stim_num))]
                    cell_entry = expt_entry.cells[stim['stimulationPoint']['name']]

                    for i in range(len(rec.pulse_start_times)):
                    # Record information about all pulses, including test pulse.
                        #t0, t1 = pulse.meta['pulse_edges']
                        #resampled = pulse['reporter'].resample(sample_rate=20000)
                        pulse_entry = db.StimPulse(
                            recording=rec_entry,
                            cell=cell_entry,
                            pulse_number=i, #pulse.meta['pulse_n'],
                            onset_time=rec.pulse_start_times[i], #t0,
                            amplitude=power_cal.convert_voltage_to_power(rec.pulse_power()[i], timestamp_to_datetime(expt_entry.acq_timestamp), expt_entry.rig_name), ## need to fill in laser/objective correctly
                            #data=resampled.data,
                            #data_start_time=resampled.t0,
                            #wavelength,
                            #light_source,
                            position=stim['stimPos'],
                            #position_offset=stim['offset'],
                            device_name=rec.device_id,
                            #qc_pass=None
                        )
                        session.add(pulse_entry)
                        pulse_entries[i] = pulse_entry

                elif 'TTL' in rec.device_id:
                    print('passing device:', rec.device_id)
                    pass
                    
                else:
                    raise Exception('need to figure out recording type')


    def job_records(self, job_ids, session):
        """Return a list of records associated with a list of job IDs.
        
        This method is used by drop_jobs to delete records for specific job IDs.
        """
        # only need to return from syncrec table; other tables will be dropped automatically.
        db = self.database
        return session.query(db.SyncRec).filter(db.SyncRec.experiment_id==db.Experiment.id).filter(db.Experiment.ext_id.in_(job_ids)).all()

    def ready_jobs(self):
        """Return an ordered dict of all jobs that are ready to be processed (all dependencies are present)
        and the dates that dependencies were created.
        """
        db = self.database
        
        # All experiments and their creation times in the DB
        expt_module = self.pipeline.get_module('opto_experiment')
        expts = expt_module.finished_jobs()
        
        # Look up nwb file locations for all experiments
        session = db.session()
        expt_recs = session.query(db.Experiment.ext_id, db.Experiment.storage_path, db.Experiment.ephys_file).filter(db.Experiment.ephys_file != None).all()
        expt_paths = {rec.ext_id: rec for rec in expt_recs}
        session.rollback()
        
        # Return the greater of NWB mod time and experiment DB record mtime
        ready = OrderedDict()
        for expt_id, (expt_mtime, expt_success) in expts.items():
            if expt_id not in expt_paths or expt_success is False:
                # no NWB file; ignore
                continue
            rec = expt_paths[expt_id]
            ephys_file = os.path.join(config.synphys_data, rec.storage_path, rec.ephys_file)
            nwb_mtime = timestamp_to_datetime(os.stat(ephys_file).st_mtime)
            ready[rec.ext_id] = max(expt_mtime, nwb_mtime)
        return ready
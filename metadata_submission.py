import os, yaml, re
from datetime import datetime, timedelta
from collections import OrderedDict
import acq4
import config
import lims
import constants


class ExperimentMetadataSubmission(object):
    """Handles final metadata QC before submission.
    
    Outputs:
    
    - pipettes.yml - list of pipettes with patch status
    - file_manifest.yml - list of files selected for LIMS submission
    - acq4 .index files updated with file categories, LIMS specimen data, rig name.
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
        
        site_info = self.site_dh.info()
        slice_dh = self.site_dh.parent()
        slice_info = slice_dh.info()
        expt_dh = slice_dh.parent()
        expt_info = expt_dh.info()
        
        spec_name = slice_info['specimen_ID'].strip()
        self.spec_info = lims.specimen_info(spec_name)

        
        # Do all categorized files actually exist?
        # Was a file categorized differently on a previous submission?
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
                if fh.isFile() and relname not in all_files and fh.ext() not in ['.pxp']:
                    warnings.append("file %s is not mentioned in metadata." % relname)
        
        # Do we already have pipette information submitted?
        self.pip_file = self.site_dh['pipettes.yml'].name()
        if os.path.exists(self.pip_file):
            yaml_pips = yaml.load(open(self.pip_file, 'rb'))
            if yaml_pips != self.pipettes:
                warnings.append('pipette metadata file will be overwritten.')
        
        # make sure we have at most one nwb file to submit
        source_nwb_files = [f['path'] for f in self.files if f['category'] == 'MIES physiology']
        if len(source_nwb_files) == 0:
            warnings.append("No NWB files specified")
        if len(source_nwb_files) > 1:
            errors.append("%d NWB files specified (should be 0 or 1)" % len(source_nwb_files))
        
        # make sure we have at most one multipatch log file to submit
        source_log_files = [f['path'] for f in self.files if f['category'] == 'Multipatch log']
        if len(source_log_files) == 0:
            warnings.append("No MultiPatch log files specified")
        if len(source_log_files) > 1:
            errors.append("%d MultiPatch log files specified (should be 0 or 1)" % len(source_log_files))

        
        # Sanity checks on experiment metadata:
        
        
        # Is rig name set and known?
        if config.rig_name not in ('MP0', 'MP1', 'MP2', 'MP3', 'MP4', 'MP5'):
            errors.append('Unrecognized rig name "%s"' % config.rig_name)
        
        # Is ACSF set?
        acsf = expt_info.get('solution', None)
        if acsf not in constants.ACSF_RECIPES:
            errors.append('Unrecognized ACSF recipe: "%s"' % acsf)
        
        # Is temperature set?
        temp = expt_info.get('temperature', '')
        m = re.match(r'(rt|\d{2}\s*c?)', temp.lower())
        if m is None:
            errors.append('Unrecognized temperature: "%s"' % temp)
        
        # Is specimen from today?
        
        
        # Sanity checks on pipette metadata:
        
        for pip_id, pip in self.pipettes.items():
            # Did we specify a known dye for each pipette?
            if pip['internal_dye'] not in constants.INTERNAL_DYES:
                errors.append('Pipette %d has unrecognized dye "%s"' % (pip_id, pip['internal_dye']))
            
            # Did we specify a known internal for each pipette?
            if pip['internal_solution'] not in constants.INTERNAL_RECIPES:
                errors.append('Pipette %d has unrecognized internal "%s"' % (pip_id, pip['internal_solution']))
            
            # Does the selected dye overlap with cre reporters?
        
        # If slice was not fixed, don't attempt LIMS submission
        try:
            sid = lims.specimen_id_from_name(spec_name)
        except ValueError as err:
            errors.append(err.message)
            sid = None

        if slice_info['plate_well_ID'] != 'not fixed' and sid is not None:
            self._check_lims(errors, warnings, spec_name, sid, site_info, slice_info)
        
        return errors, warnings
        
    def _check_lims(self, errors, warnings, spec_name, sid, site_info, slice_info):
        cw_name = self.spec_info['carousel_well_name']
        
        # LIMS upload will fail if the specimen has not been given an ephys roi plan.
        roi_plans = lims.specimen_ephys_roi_plans(spec_name)
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
        
        if self.spec_info['flipped'] not in (True, False):
            errors.append("Must set flipped field for this specimen: %s" % lims_edit_href)
        # histology well name should look something like "multi_170911_21_A01"
        hist_well = self.spec_info['histology_well_name']
        if hist_well is None:
            errors.append("Missing histology well name in LIMS: %s" % lims_edit_href)        
        else:
            m = re.match(r'multi_(\d{2})(\d{2})(\d{2})_(2[1-9])_([A-C]0[1-4])', hist_well)
            if m is None:
                errors.append("Histology well name appears to be incorrectly formatted: %s" % lims_edit_href)        
            else:
                yy, mm, dd, plate_n, well = m.groups()[:5]
                plate_date = datetime(2000+int(yy), int(mm), int(dd))
                site_date = datetime.fromtimestamp(site_info['__timestamp__'])
                # find the most recent Monday
                expected_plate_date = site_date - timedelta(days=site_date.weekday())
                if abs((expected_plate_date - plate_date).total_seconds()) > 7*24*3600:
                    # error if more than a week out of sync
                    errors.append("Histology well date is %s%s%s; expected %s: %s" % (yy, mm, dd, expected_plate_date.strftime('%y%m%d'), lims_edit_href))
                if expected_plate_date.date() != plate_date.date():
                    # warning if the date is not exactly as expected
                    warnings.append("Histology well date is %s%s%s; expected %s: %s" % (yy, mm, dd, expected_plate_date.strftime('%y%m%d'), lims_edit_href))
                    
                if int(plate_n) > 24:
                    warnings.append("Histology plate number %s is probably too high. %s" % (plate_n, lims_edit_href))
                elif int(plate_n) > 22:
                    warnings.append("Histology plate number %s might be too high? %s" % (plate_n, lims_edit_href))

        # Check carousel ID matches the one in LIMS
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

    def summary(self):
        summ = OrderedDict()
        summ['file categories'] = self.files
        summ['pipettes'] = self.pipettes
        
        return {'metadata': summ}
    
    def submit(self):
        slice_dh = self.site_dh.parent()
        expt_dh = slice_dh.parent()
        
        # Set rig name from config if it has not already been written
        if 'rig_name' not in expt_dh.info():
            expt_dh.setInfo(rig_name=config.rig_name)
            
        # Write LIMS info to top-level directory handle.
        # This is just for convenience when browsing the raw data.
        expt_dh.setInfo(LIMS_specimen_info=self.spec_info)
            
        for file_info in self.files:
            if file_info['category'] == 'ignore':
                continue
            fh = slice_dh[file_info['path']]
            fh.setInfo(category=file_info['category'])
        
        # Generate yml file with pipette/cell information. Operators can edit this
        # information later on..
        yaml.dump(self.pipettes, open(self.pip_file, 'wb'), default_flow_style=False, indent=4)

        # Write an manifest file describing the data that should be uploaded to LIMS
        manifest_file = self.site_dh['file_manifest.yml'].name()
        yaml.dump(self.files, open(manifest_file, 'wb'), default_flow_style=False, indent=4)



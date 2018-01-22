import os, yaml, re
from datetime import datetime, timedelta
from collections import OrderedDict
import acq4

from . import config, lims, constants, genotypes, experiment_list
from . import yaml_local  # adds support for OrderedDict


class ExperimentMetadataSubmission(object):
    """Handles final metadata QC before submission.
    
    Outputs:
    
    - pipettes.yml - list of pipettes with patch status
    - file_manifest.yml - list of files selected for LIMS submission
    - acq4 .index files updated with file categories, LIMS specimen data, rig name.

    If possible, some fields in pipettes.yml are filled from old-format metadata files
    """
    message = "Copying submission metadata back to index files"

    _expt_list = None  # used for importing old metadata
    
    def __init__(self, site_dh, files, pipettes):
        self.site_dh = site_dh
        self.files = files
        self.pipettes = pipettes

    @property
    def old_expt_list(self):
        """Singleton ExperimentList used to access old metadata information for import
        """
        if self._expt_list is None:
            ExperimentMetadataSubmission._expt_list = experiment_list.ExperimentList(cache='expts_cache.pkl')
        return self._expt_list

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
        # Are file sizes appropriate?
        max_sizes = {
            'MIES physiology': 5000,
            'Multipatch log': 1,
            'slice anatomy': 10,
            'slice quality stack': 5000,
            'recording site': 5000,
        }
        for file_info in self.files:
            fh = slice_dh[file_info['path']]
            if not fh.isFile():
                errors.append("file %s does not exist" % fh.name())
                continue
            info = fh.info()
            if 'category' in info and info['category'] != file_info['category']:
                warnings.append("category for file %s will change %s => %s" % 
                                (info['category'], file_info['category']))
            
            max_size = max_sizes.get(file_info['category'], 0)
            size = int(os.stat(fh.name()).st_size / 1e6)
            if max_size > 0 and size > max_size:
                errors.append('File "%s" is too large (%dMB > %dMB max for category "%s")' % 
                              (fh.name(), size, max_size, file_info['category']))
            
        
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
        have_nwb = len(source_nwb_files) > 0
        if len(source_nwb_files) == 0:
            warnings.append("No NWB files specified")
        if len(source_nwb_files) > 1:
            errors.append("%d NWB files specified (should be 0 or 1)" % len(source_nwb_files))
        
        # make sure we have one multipatch log file to submit if an nwb is present
        # or 0-1 log files if there is no nwb
        source_log_files = [f['path'] for f in self.files if f['category'] == 'Multipatch log']
        if have_nwb and len(source_log_files) == 0:
            # Occasionally user forgets to record this data; usually we can work with it
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
        

        # Attempt import of old-format metadata
        self._import_old_metadata(site_info, warnings, errors)

        return errors, warnings
        
    def _check_lims(self, errors, warnings, spec_name, sid, site_info, slice_info):
        cw_name = self.spec_info['carousel_well_name']
        
        # LIMS upload will fail if the specimen has not been given an ephys roi plan.
        accepted_roi_plans = {
            'mouse': 'Synaptic Physiology ROI Plan', 
            'human': 'Synaptic Physiology Human ROI Plan'
        }[self.spec_info['organism']]
        roi_plans = lims.specimen_ephys_roi_plans(spec_name)
        lims_edit_href = '<a href="http://lims2/specimens/{sid}/edit">http://lims2/specimens/{sid}/edit</a>'.format(sid=sid)
        if len(roi_plans) == 0:
            errors.append('Specimen has no ephys roi plan. Edit:' + lims_edit_href)
        elif len(roi_plans) == 1 and roi_plans[0]['name'] not in accepted_roi_plans:
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

        # check specimen age
        if self.spec_info['age'] is None:
            warnings.append("Donor age is nto set in LIMS.")

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
            
        # Write LIMS info to top-level and slice directory handles.
        # This is just for convenience when browsing the raw data.
        slice_info = self.spec_info.copy()
        donor_keys = ['organism', 'age', 'date_of_birth', 'genotype', 'weight', 'sex']
        donor_info = {k:slice_info.pop(k) for k in donor_keys}
        expt_dh.setInfo(LIMS_donor_info=donor_info)
        slice_dh.setInfo(LIMS_specimen_info=slice_info)
            
        # Write category labels for each image
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

    def _import_old_metadata(self, site_info, warnings, errors):
        expts = self.old_expt_list
        uid = '%0.02f' % site_info['__timestamp__']
        try:
            genotype = genotypes.Genotype(self.spec_info['genotype'])
        except Exception:
            genotype = None
            warnings.append('Error parsing genotype "%s"; will not be able to import old cre label calls' % (self.spec_info['genotype']))
        
        try:
            expt = expts[uid]
        except KeyError:
            warnings.append("Could not import any old-format metadata (uid=%s)." % uid)
            return

        connections = expt.connections[:]
        for cell in expt.cells.values():
            if cell.cell_id not in self.pipettes:
                warnings.append("Old-format metadata contains a cell ID %s, but this does not exist in the current submission." % cell.cell_id)
                continue

            pid = cell.cell_id
            pip = self.pipettes[pid]

            # import labels
            labels = {}
            for label,pos in cell._raw_labels.items():
                if pos not in ['', '+', '-', '+?', '-?', 'x']:
                    warnings.append('Pipette %d: ignoring old label "%s" because the value "%s" is unrecognized' % (pid, label, pos))
                    continue

                # biocytin or fluorophore
                if label == 'biocytin':
                    labels[label] = pos
                elif label in constants.FLUOROPHORES:
                    color = constants.FLUOROPHORES[label]
                    labels[color] = pos
                # human_L layer call
                elif label.startswith('human_L'):
                    if pos == 'x':
                        continue
                    layer = label[7:]
                    if pip['target_layer'] == '':
                        pip['target_layer'] = layer
                        warnings.append("Pipette %d: imported layer %s from old metadata." % (pid, layer))
                    elif pip['target_layer'] != layer:
                        warnings.append("Pipette %d: old metadata layer %s conflicts with current layer: %s." % (pid, layer, pip['target_layer']))
                # cre type; convert to color(s)
                elif label in constants.ALL_CRE_TYPES:
                    if genotype is None:
                        warnings.append("Pipette %d: ignoring old cre label %s" % (pid, label))
                        continue
                    if label not in genotype.drivers():
                        warnings.append("Pipette %d: old cre label %s is not in genotype!" % (pid, label))
                        continue
                    for color in genotype.colors(driver=label):
                        if color in labels:
                            warnings.append("Pipette %d: color %s is specified twice!" % (pid, color))
                        labels[color] = pos
                else:
                    warnings.append("Pipette %d: old metadata has unrecognized label: %s." % (pid, label))
            
            # now make sure there are no conflicts
            for label, pos in labels.items():
                val = pip['cell_labels'].get(label, '')
                if val == '':
                    pip['cell_labels'][label] = pos
                    warnings.append("Pipette %d: imported label %s=%s from old metadata." % (pid, label, pos))
                elif val != pos:
                    warnings.append("Pipette %d: old metadata laybel %s=%s conflicts with current value: %s." % (pid, label, pos, val))

            # import old QC (we can be more lax here because this is deprecated data)
            pip['cell_qc'] = {'holding': cell.holding_qc, 'access': cell.access_qc, 'spiking': cell.spiking_qc}

            # import connections
            cell_connects = []
            for pre, post in connections[:]:
                if pre != pid:
                    continue
                connections.remove((pre, post))
                cell_connects.append(post)
            if pip['synapse_to'] is None:
                pip['synapse_to'] = cell_connects
                warnings.append("Pipette %d: imported connections %s from old metadata." % (pid, cell_connects))
            else:
                if list(sorted(pip['synapse_to'])) != list(sorted(cell_connects)):
                    warnings.append("Pipette %d: old metadata connections %s conflicts with current value: %s." % (pid, cell_connects, pip['synapse_to']))

        if len(connections) > 0:
            warnings.append("Could not import old metadata connections:" % (connections))

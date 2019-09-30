import os, sys, glob, json
import yaml
from collections import OrderedDict
from acq4.util.DataManager import getDirHandle
from aisynphys import config, lims


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
        self.spec_id = sid
        
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

    def status(self):
        ts = self.meta['acq_timestamp']
        return lims.expt_submissions(self.spec_id, ts)


def find_submittable_expts():
    """Search synphys data storage for experiments that are ready to be submitted to LIMS.
    """
    
    # Start by finding all site paths that have an nwb file and a file_manifest.yml
    all_nwbs = glob.glob(os.path.join(config.synphys_data, '*', 'slice_*', 'site_*', '*.nwb'))
    all_sites = OrderedDict()    
    for nwb in all_nwbs:
        path = os.path.dirname(nwb)
        if 'file_manifest.yml' in os.listdir(path):
            all_sites[path] = 1
    
    # filter out anything that has been submitted already
    ready_sites = []
    for path in all_sites.keys():
        site_dh = getDirHandle(path)
        acq_ts = site_dh.info()['__timestamp__']
        spec_name = site_dh.parent().info()['specimen_ID'].strip()
        spec_id = lims.specimen_id_from_name(spec_name)
        subs = lims.expt_submissions(spec_id, acq_ts)
        if len(subs) == 0:
            ready_sites.append(site_dh)

    return ready_sites
    

def submit_site(dh):
    files = yaml.load(open(dh['file_manifest.yml'].name(), 'rb'))
    
    # file manifests are written by windows, so we have to re-normalize the file names.
    for f in files:
        f['path'] = os.path.join(*f['path'].split('\\'))
        
    sub = LIMSSubmission(dh, files)
    err, warn = sub.check()
    if len(err) > 0:
        raise Exception('LIMS submission errors:' + '\n'.join(err))
    
    sub.submit()
    
    return sub
    
    
if __name__ == '__main__':
    print("Searching for sites to submit..")
    sites = find_submittable_expts()
    
    subs = []
    for i,site in enumerate(sites):
        print("=== Submit site %d/%d : %s" % (i, len(sites), site))
        try:
            sub = submit_site(site)
            subs.append(sub)
        except Exception as exc:
            sys.excepthook(*sys.exc_info())
            subs.append(exc)
            
    
    
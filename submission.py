import datetime
from database import Session, Slice
from allensdk_internal.core import lims_utilities as lims

import config

class SliceSubmission(object):
    def __init__(self, dh):
        self.dh = dh
        
        self._fields = None

    @property
    def fields(self):
        info = dh.info()
        
        # pull some metadata from LIMS
        sid = info['specimen_ID']
        query = """
            select 
              organisms.name as organism, 
              ages.days as age, 
              donors.full_genotype as genotype,
              tissue_processings.section_thickness_um as thickness,
              plane_of_sections.name as orientation
            from specimens
              left join donors on specimens.donor_id=donors.id 
              left join organisms on donors.organism_id=organisms.id
              left join ages on donors.age_id=ages.id
              left join tissue_processings on specimens.tissue_processing_id=tissue_processings.id
              left join plane_of_sections on tissue_processings.plane_of_section_id=plane_of_sections.id
            where specimens.name='%s';
        """ % sid
        r = lims.query(query)
        if len(r) != 1:
            raise Exception("LIMS lookup for specimen %s returned %d results (expected 1)" % (sid, len(r)))
        limsdata = r[0]
        
        self._fields = {
            'acq_timestamp', datetime.fromtimestamp(info['__timestamp__']),
            'species': limsdata['organism'],
            'age': limsdata['age'],
            'genotype': limsdata['genotype'],
            'orientation' limsdata['orientation'],
            'surface': info.get('surface'),
            'quality': info.get('slice quality'),
            'slice_time': None,
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
        if fields['surface'] not in ['medial', 'lateral']:
            messages.append("Warning: slice surface '%s' should have been 'medial' or 'lateral'" % surface)
        return errors, warnings
        
    def create(self):
        if len(self.check()[0]) > 0:
            raise Exception("Submission has errors; see SliceSubmission.check()")
        data = self.fields
        sl = Slice(lims_specimen_name=data['specimen_name'], surface=data['surface'],
                original_path=data['original_path'], acq_timestamp=data['acq_timestamp'],
                submission_data=data['image_files'], quality=data['slice_quality'])
        return sl
        


def submit_slice(data):
    """Submit information about a new slice to LIMS and the internal analysis DB.
    
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
    session = Session()
    session.commit()
    
    

def submit_experiment(data):
    """Submit a new experiment to LIMS and the internal analysis DB.
    
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
    
    

from database import Session, Slice


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
    sl = Slice(lims_specimen_name=data['specimen_name'], surface=data['surface'],
               original_path=data['original_path'], acq_timestamp=data['acq_timestamp'],
               submission_data=data['image_files'], quality=data['slice_quality'])
    session.add(sl)
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
    
    

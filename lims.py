import re
from allensdk_internal.core import lims_utilities as lims


def specimen_info(specimen_name):
    """Return a dictionary of information about a specimen queried from LIMS.
    
    Also generates information about the hemisphere and which side of the slice
    was patched.
    
    Returns
    -------
    organism : "mouse" or "human"
    genotype : the full genotype of the donor as recorded in labtracks
    age : age of specimen in days at time of sectioning
    date_of_birth : donor's date of birth
    plane_of_section : 'coronal' or 'sagittal'
    hemisphere : 'left' or 'right'
    thickness : speimen slice thickness (unscaled meters)
    sectioning_mount_side : the side of the tissue that was mounted during
        sectioning
    flipped : boolean; if True, then the slice was flipped relative to its 
        blockface image during recording
    exposed_surface : The surface that was exposed during the experiment (right, 
        left, anterior, or posterior)
    section_instructions : description of the slice angle and target region
        used for sectioning
    section_number : indicates the order this slice was sectioned (1=first)
    """
    
    # Query all interesting information about this specimen from LIMS
    sid = specimen_name.strip()
    query = """
        select 
            organisms.name as organism, 
            ages.days as age,
            donors.date_of_birth as date_of_birth,
            donors.full_genotype as genotype,
            tissue_processings.section_thickness_um as thickness,
            tissue_processings.instructions as section_instructions,
            plane_of_sections.name as plane_of_section,
            flipped_specimens.name as flipped
        from specimens
            left join donors on specimens.donor_id=donors.id 
            left join organisms on donors.organism_id=organisms.id
            left join ages on donors.age_id=ages.id
            left join tissue_processings on specimens.tissue_processing_id=tissue_processings.id
            left join plane_of_sections on tissue_processings.plane_of_section_id=plane_of_sections.id
            left join flipped_specimens on flipped_specimens.id = specimens.flipped_specimen_id
        where specimens.name='%s';
    """ % sid
    r = lims.query(query)
    if len(r) != 1:
        raise Exception("LIMS lookup for specimen '%s' returned %d results (expected 1)" % (sid, len(r)))
    rec = r[0]
    
    # convert thickness to unscaled
    rec['thickness'] = rec['thickness'] * 1e-6
    # convert organism to more easily searchable form
    rec['organism'] = {'Mus musculus': 'mouse', 'Homo Sapiens': 'human'}[rec['organism']]
    # convert flipped to bool
    rec['flipped'] = {'flipped': True, 'not flipped': False, 'not checked': None}[rec['flipped']]
    
    # Parse the specimen name to extract more information about the plane of section.
    # Format is:  
    #    driver1;reporter1;driver2;reporter2-AAAAAA.BB.CC
    #        AAAAAA = donor ID
    #            BB = slice number
    #            CC = orientation and hemisphere
    m = re.match(r'(.*)-(\d{6,7})(\.(\d{2}))(\.(\d{2}))$', sid)
    if m is None:
        raise Exception('Could not parse specimen name: "%s"' % sid)
    
    rec['section_number'] = int(m.groups()[3])
    
    # The last number contains information about the orientation and hemisphere
    orientation_num = m.groups()[5]
    plane, hem, mount = {
        '01': ('coronal', 'left', 'anterior'),
        '02': ('coronal', 'right', 'anterior'),
        '03': ('coronal', 'left', 'posterior'),
        '04': ('coronal', 'right', 'posterior'),
        '05': ('sagittal', 'left', 'right'),
        '06': ('sagittal', 'right', 'left'),
    }[orientation_num]
    
    if plane != rec['plane_of_section']:
        raise ValueError("Slice orientation from specimen name (.%s=%s) does not match plane_of_section recorded in LIMS (%s)" % (orientation_num, plane, rec['plane_of_section']))
    
    rec['hemisphere'] = hem
    rec['sectioning_mount_side'] = mount
    
    # decide which surface was patched
    if rec['flipped'] is True:
        rec['exposed_surface'] = mount
    elif rec['flipped'] is False:
        rec['exposed_surface'] = {'anterior': 'posterior', 'posterior': 'anterior', 'right': 'left', 'left': 'right'}[mount]
    else:    
        rec['exposed_surface'] = None
    
    return rec
    
    
def specimen_images(specimen_name):
    """Return a list of (image ID, treatment) pairs for a specimen.
    
    """
    q = """
        select sub_images.id, treatments.name from specimens 
        join image_series on image_series.specimen_id=specimens.id 
        join sub_images on sub_images.image_series_id=image_series.id
        join images on images.id = sub_images.image_id
        left join treatments on treatments.id = images.treatment_id
        where specimens.name='%s';
        """ % specimen_name
    r = lims.query(q)
    return [(rec['id'], rec['name']) for rec in r]

        
def submit_expt(spec_id, nwb_file, json_file):
    import limstk.LIMStk as limstk
    #limstk.init_log()

    """
    This is an example of creating a trigger file and copying files for multi-patch.
    lims_scheduler_d.py should be running to receive requests from this library.

    When creating the session, HWBIgor and metadata are values you would be passing in (per confluence spec) rather than
    retrieving from lims2/ but 'id' will be the value you are passing to lims2/ to get the trigger directory and roi plans
    back from.

    These are defined in the lims_config.yml
    """

    # notice the "' notation when adding the filename - this is to allow the ' to appear in the trigger file
    # 'multipatch' is a lookup key for your specification but,
    # id, NWBIgor and metadata are key-value pairs you want to see in the trigger file but that do not come from lims
    # they are also defined in the limstk_config.
    lims_session = limstk.Session(
        'multipatch',  # defined in the limstk_config file
        id=spec_id,
        NWBIgor="'/allen/programs/celltypes/production/incoming/mousecelltypes/multipatch-test.nwb'",
        metadata="'/allen/programs/celltypes/production/incoming/mousecelltypes/multipatch-test.json'")

    # Because there could be multiple plans returned on this query, the user has to determine which plan id is correct
    # For these situations, there are 'manual' requests like specimens_by_id and specimens_by_well_name
    resp = lims_session.request('specimens_by_id', id=spec_id)

    # This data can be manually added to the trigger data
    lims_session.trigger_data['id'] = resp['ephys_specimen_roi_plans'][0]['id']

    # enumerate the files you'd like to copy over
    lims_session.add_to_manifest(json_file, dst_filename='multipatch-test.json')
    lims_session.add_to_manifest(nwb_file, dst_filename='multipatch-test.nwb')

    # you could optionally copy a file over with a new name
    # lims_session.add_to_manifest('c:/myFile', dst_filename = 'newFilename') <--- no path necessary on dst_filename

    # to finish out, schedule the session
    lims_session.commit_manifest(trigger_file='%d.mp' % spec_id)


if __name__ == '__main__':
    spec_name = "Ntsr1-Cre_GN220;Ai14-349905.03.06"
    recs = lims.query("select id from specimens where name='%s'" % spec_name)
    spec_id = recs[0]['id']
    print (spec_id)
    submit_expt(spec_id, 'lims_test.nwb', 'lims_test.json')
    
    
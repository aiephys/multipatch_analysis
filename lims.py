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
    sid = specimen_name
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
        raise Exception("LIMS lookup for specimen %s returned %d results (expected 1)" % (sid, len(r)))
    rec = r[0]
    
    # convert thickness to unscaled
    rec['thickness'] = rec['thickness'] * 1e-6
    # convert organism to more easily searchable form
    rec['organism'] = {'Mus musculus': 'mouse', 'Homo Sapiens': 'human'}[rec['organism']]
    # convert flipped to bool
    rec['flipped'] = {'flipped': True, 'not flipped': False, 'not checked': None}[rec['flipped']]
    
    # Parse the specimen name to extract more information about the plane of section
    m = re.match(r'(.*)-(\d{6,7})(\.(\d{2}))(\.(\d{2}))$', sid)
    if m is None:
        raise Exception('Could not parse specimen name: "%s"' % sid)
    
    rec['section_number'] = int(m.groups()[3])
    
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

        
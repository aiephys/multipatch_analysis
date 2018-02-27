from __future__ import print_function
import os, re, json
from allensdk_internal.core import lims_utilities as lims


def specimen_info(specimen_name=None, specimen_id=None):
    """Return a dictionary of information about a specimen queried from LIMS.
    
    Also generates information about the hemisphere and which side of the slice
    was patched.
    
    Returns
    -------
    organism : "mouse" or "human"
    age : age of specimen in days at time of sectioning
    date_of_birth : donor's date of birth
    genotype : the full genotype of the donor as recorded in labtracks
    weight : weight in grams
    sex : 'M', 'F', or 
    plane_of_section : 'coronal' or 'sagittal'
    hemisphere : 'left' or 'right'
    thickness : speimen slice thickness (unscaled meters)
    section_instructions : description of the slice angle and target region
        used for sectioning
    flipped : boolean; if True, then the slice was flipped relative to its 
        blockface image during recording
    sectioning_mount_side : the side of the tissue that was mounted during
        sectioning
    exposed_surface : The surface that was exposed during the experiment (right, 
        left, anterior, or posterior)
    section_number : indicates the order this slice was sectioned (1=first)
    """
    
    # Query all interesting information about this specimen from LIMS
    query = """
        select 
            organisms.name as organism, 
            ages.days as age,
            donors.date_of_birth as date_of_birth,
            donors.full_genotype as genotype,
            donors.weight as weight,
            genders.name as sex,
            tissue_processings.section_thickness_um as thickness,
            tissue_processings.instructions as section_instructions,
            plane_of_sections.name as plane_of_section,
            flipped_specimens.name as flipped,
            specimens.histology_well_name as histology_well_name,
            specimens.carousel_well_name as carousel_well_name,
            specimens.parent_id as parent_id,
            specimens.name as specimen_name,
            specimens.id as specimen_id
        from specimens
            left join donors on specimens.donor_id=donors.id 
            left join organisms on donors.organism_id=organisms.id
            left join ages on donors.age_id=ages.id
            left join genders on donors.gender_id=genders.id
            left join tissue_processings on specimens.tissue_processing_id=tissue_processings.id
            left join plane_of_sections on tissue_processings.plane_of_section_id=plane_of_sections.id
            left join flipped_specimens on flipped_specimens.id = specimens.flipped_specimen_id
    """
    if specimen_name is not None:
        sid = specimen_name.strip()
        query += "where specimens.name='%s';" % sid
    elif specimen_id is not None:
        sid = specimen_id
        query += "where specimens.id='%d';" % sid
    else:
        raise ValueError("Must specify specimen name or ID")
        
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
    # Mouse format is:  
    #    driver1;reporter1;driver2;reporter2-AAAAAA.BB.CC
    #        AAAAAA = donor ID
    #            BB = slice number
    #            CC = orientation and hemisphere
    spec_name = rec['specimen_name']
    if rec['organism'] == 'mouse':
        m = re.match(r'(.*)(-(\d{6,7}))?(\.(\d{2}))(\.(\d{2}))$', spec_name)
        if m is None:
            raise Exception('Could not parse mouse specimen name: "%s"' % spec_name)
        
        rec['section_number'] = int(m.groups()[4])
        
        # The last number contains information about the orientation and hemisphere
        orientation_num = m.groups()[6]
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
            
    # Human format is:
    #   Haa.bb.ccc.dd.ee.ff
    elif rec['organism'] == 'human':
        m = re.match(r'H(\d+)\.(\d+)\.(\d+)\.(\d+)\.(\d+)(\.(\d+))?$', spec_name)
        if m is None:
            raise Exception('Could not parse human specimen name: "%s"' % spec_name)
        rec['hemisphere'] = None
        rec['sectioning_mount_side'] = None
        rec['exposed_surface'] = None
        rec['human_donor_site'] = int(m.groups()[1])
        rec['human_donor_number'] = int(m.groups()[2])
        rec['block_number'] = int(m.groups()[3])
        rec['section_number'] = int(m.groups()[4])
        rec['subsection_number'] = None if m.groups()[6] is None else int(m.groups()[6])
        
        
    else:
        raise Exception('Unsupported organism: "%s"' % rec['organism'])
    
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


def specimen_id_from_name(spec_name):
    """Return the LIMS ID of a specimen give its name.
    """
    recs = lims.query("select id from specimens where name='%s'" % spec_name)
    if len(recs) == 0:
        raise ValueError('No LIMS specimen named "%s"' % spec_name)
    return recs[0]['id']


def specimen_ephys_roi_plans(spec_name):
    """Return a list of all ephys roi plans for this specimen.
    """
    sid = specimen_id_from_name(spec_name)
    recs = lims.query("""
        select 
            ephys_roi_plans.id as ephys_roi_plan_id,
            ephys_specimen_roi_plans.id as ephys_specimen_roi_plan_id,
            ephys_roi_plans.name as name
        from 
            ephys_specimen_roi_plans
            join ephys_roi_plans on ephys_specimen_roi_plans.ephys_roi_plan_id=ephys_roi_plans.id
        where 
            ephys_specimen_roi_plans.specimen_id=%d
    """ % sid)
    return recs


def cell_cluster_ids(spec_id):
    recs = lims.query("select id from specimens where specimens.parent_id=%d" % spec_id)
    return [rec['id'] for rec in recs]


def cell_cluster_data_paths(cluster_id):
    recs = lims.query("""
        select ephys_roi_results.storage_directory 
        from specimens 
        join ephys_roi_results on ephys_roi_results.id=specimens.ephys_roi_result_id
        where specimens.id=%d
    """ % cluster_id)
    return [r['storage_directory'] for r in recs]


def specimen_metadata(spec_id):
    recs = lims.query("select data from specimen_metadata where specimen_id=%d" % spec_id)
    if len(recs) == 0:
        return None
    meta = recs[0]['data']
    if meta == '':
        return None
    return json.loads(meta)  # unserialization corrects for a LIMS bug; we can remove this later.


def filename_base(specimen_id, acq_timestamp):
    """Return a base filename string to be used for all LIMS uploads.
    """
    return "synphys-%d-%s" % (specimen_id, acq_timestamp)


def submit_expt(spec_name, acq_timestamp, nwb_file, json_file):
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
    spec_id = specimen_id_from_name(spec_name)

    # notice the "' notation when adding the filename - this is to allow the ' to appear in the trigger file
    # 'multipatch' is a lookup key for your specification but,
    # id, NWBIgor and metadata are key-value pairs you want to see in the trigger file but that do not come from lims
    
    filebase = filename_base(spec_id, acq_timestamp)
    
    incoming_files = dict(
        NWBIgor="/allen/programs/celltypes/production/incoming/mousecelltypes/%s.nwb" % filebase,
        metadata="/allen/programs/celltypes/production/incoming/mousecelltypes/%s.json" % filebase,
    )
    
    # they are also defined in the limstk_config.
    lims_session = limstk.Session(
        'multipatch',  # defined in the limstk_config file
        id=spec_id,
        NWBIgor="'" + incoming_files['NWBIgor'] + "'",
        metadata="'" + incoming_files['metadata'] + "'"
    )

    # Because there could be multiple plans returned on this query, the user has to determine which plan id is correct
    # For these situations, there are 'manual' requests like specimens_by_id and specimens_by_well_name
    resp = lims_session.request('specimens_by_id', id=spec_id)

    # This data can be manually added to the trigger data
    lims_session.trigger_data['id'] = resp['ephys_specimen_roi_plans'][0]['id']

    # enumerate the files you'd like to copy over
    lims_session.add_to_manifest(json_file, dst_filename='%s.json' % filebase)
    lims_session.add_to_manifest(nwb_file, dst_filename='%s.nwb' % filebase)

    # to finish out, schedule the session
    lims_session.commit_manifest(trigger_file='%s.mp' % filebase)

    return incoming_files


def expt_submissions(spec_id, acq_timestamp):
    """Return information about the status of each submission found for an
    experiment, identified by its specimen ID and experiment uid.
    """
    submissions = []
    filebase = filename_base(spec_id, acq_timestamp)
    
    # Do we have incoming files?
    incoming_path = '/allen/programs/celltypes/production/incoming/mousecelltypes'
    incoming_nwb = os.path.join(incoming_path, filebase + '.nwb')
    incoming_trigger = os.path.join(incoming_path, 'trigger', '%s.mp' % filebase)
    failed_trigger = os.path.join(incoming_path, 'failed_trigger', '%s.mp' % filebase)
    
    if os.path.exists(incoming_nwb):
        if os.path.exists(incoming_trigger):
            submissions.append(("trigger pending", incoming_trigger))
        if os.path.exists(failed_trigger):
            error = open(failed_trigger+'.err', 'r').read()
            submissions.append(("trigger failed", failed_trigger, error))
            
    # Anything in LIMS already?
    cluster_ids = cell_cluster_ids(spec_id)
    for cid in cluster_ids:
        meta = specimen_metadata(cid)
        if meta is not None and meta['acq_timestamp'] == acq_timestamp:
            data_path = cell_cluster_data_paths(cid)
            submissions.append(("succeeded", cid, data_path))
    
    return submissions
    

if __name__ == '__main__':
    # testing specimen
    spec_name = "Ntsr1-Cre_GN220;Ai14-349905.03.06"
    spec_id = specimen_id_from_name(spec_name)
    cluster_ids = cell_cluster_ids(spec_id)
    print("Slice:", spec_id)
    print(specimen_info(spec_name))
    for cid in cluster_ids:
        print("  %d %s" % (cid, specimen_metadata(cid)))
    #print(cell_cluster_data_paths(cluster_ids[0]))
    print("")
    for sub in expt_submissions(spec_id, 1505768693.087):
        print(sub)

    #submit_expt(spec_id, 'lims_test.nwb', 'lims_test.json')
    
    
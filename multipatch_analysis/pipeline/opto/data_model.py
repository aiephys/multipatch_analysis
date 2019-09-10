from multipatch_analysis import lims
import re


def find_lims_specimen_ids(slice_dh):
    """Return a list of lims specimen_ids that match the metainfo in the day and slice .index files.
    Search order:
        1) slice.specimen_ID
        2) day.animal_ID + slice.slice_number

    """

    info = slice_dh.info()
    parent_info = slice_dh.parent().info()

    sid = info.get('specimen_ID', '').strip()
    if sid == '':
        slice_id = info.get('slice_number', '').strip()
        if len(slice_id) > 2:
            sid = slice_id
        else:
            animal_id = parent_info.get('animal_id', '').strip()
            if len(animal_id) == 0:
                animal_id = parent_info.get('animal_ID', '').strip()
                if len(animal_id) == 0:
                    return []
            if len(slice_id) == 1:
                slice_id = '0'+slice_id
            sid = animal_id + '.' + slice_id

    #print('sid:', sid)
    ids = lims.find_specimen_ids_matching_name(sid)
    if len(ids) == 1:
        return ids

    possible_ids = []
    for n in ids:
        r = lims.query("select specimens.name as specimen_name from specimens where specimens.id=%d"%n)
        if len(r) != 1:
            raise Exception("LIMS lookup for specimen '%s' returned %d results (expected 1)" % (str(n), len(r)))
        rec = dict(r[0])
        m = re.match(r'(.*)(-(\d{6,7}))?(\.(\d{2}))(\.(\d{2}))$', rec['specimen_name'])
        if m is not None:
            possible_ids.append(n)

    return possible_ids
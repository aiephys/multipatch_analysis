import os, sys, glob, datetime, csv, argparse
from acq4.util.DataManager import getDirHandle
from aisynphys import config
from aisynphys.constants import INTERNAL_RECIPES
from aisynphys.experiment import Experiment

slice_csv_file = 'Slice_Preparation_Record.csv'
osm_csv_file = 'MultiPatch Throughput.xlsx - Solution Making.csv'

# read all dissection times
recs = list(csv.reader(open(slice_csv_file).readlines()))
diss_columns = {
    'id_col': recs[0].index('LabTracks ID'),
    'date_col': recs[0].index('Date Sliced'),
    'sac_col': recs[0].index('Time Mouse Sacrificed'),
    'time_col': recs[0].index('heated ACSF1 incub. Comp.')}
diss_times = {rec[diss_columns['id_col']]: rec for rec in recs[1:]}

# read solution making
osm_recs = list(csv.reader(open(osm_csv_file).readlines()))
osm_columns = {
    'osm_date_col': osm_recs[0].index('Date'),
    'recipe_col': osm_recs[0].index('Recipe'),
    'osm_col': osm_recs[0].index('osmolarity')}
osm_dates = {osm_rec[osm_columns['osm_date_col']]: osm_rec for osm_rec in osm_recs[1:]}

def solution_check(dh, expt_date, columns, set_data=False):
    solution = dh.info().get('solution', None)
    osm = dh.info().get('solution_osm', None)
    date = expt_date.strftime('%#m/%#d/%Y')
    if expt_date < datetime.date(2016, 10, 18):
        return None
    try:
        osm_entry = osm_dates[date]
    except KeyError:
        print_msg = ("\tSolution date not in Multipatch_throughput sheet, check notes")
        return print_msg

    recipe = osm_entry[columns['recipe_col']]
    sheet_osm = osm_entry[columns['osm_col']]
    if solution is None:
        if set_data is True:
            dh.setInfo(solution=recipe)
            return None
        else:
            print_msg = ("\tSet aCSF: %s" % recipe)
            return print_msg

    calcium = solution.split('m')[0]
    if calcium not in recipe:
        print_msg = ("\taCSF disagreement: %s != %s" % (solution, recipe))
    else:
        print_msg = ''

    if osm is None:
        if set_data is True:
            dh.setInfo(solution_osm=sheet_osm)
            return None
        else:
            print_msg += ("\tSet osmolarity: %s" % sheet_osm)
            return print_msg

    if sheet_osm != osm:
        print_msg += ("\tOsmolarity disagreement: %s != %s" % (osm, sheet_osm))
        return print_msg
    return None


def dissection_check(dh, sub_id, expt_date, columns, set_data=False):
    if sub_id is None:
        print_msg = ("\tNo animal_ID for %s" % dh.path)
        return print_msg
    if sub_id.startswith('H'):
        return None
    dis_rec = diss_times.get(sub_id, None)
    if dis_rec is None:
        print_msg = ("\tNo tissue processing record")
        return print_msg

    date = dis_rec[columns['date_col']]
    sac = dis_rec[columns['sac_col']]
    time = dis_rec[columns['time_col']]
    sac_time = datetime.datetime.strptime(date + ' ' + sac, '%m/%d/%Y %I:%M %p')
    if time in (None, 'N/A', ''):
        print_msg = ("\tNo time record, sacrifice time was %s\n" % sac)
        return print_msg
    if time[-2:] in ('AM', 'PM'):
        dis_time = datetime.datetime.strptime(date + ' ' + time, '%m/%d/%Y %I:%M %p')
    else:
        if date.endswith( 'N/A'): date = date[:-4]
        dis_time = datetime.datetime.strptime(date + ' ' + time, '%m/%d/%Y %I:%M')
        if dis_time < sac_time:
            dis_time = dis_time + datetime.timedelta(seconds=12 * 3600)

    if (dis_time - sac_time).total_seconds() > 3600:
        print_msg = ("\tSac. / dissection times for are too far apart")
        return print_msg
    if dis_time <= sac_time:
        print_msg = ("\tDissection before sac")
        return print_msg

    if expt_date != dis_time.date():
        if (expt_date - dis_time.date()).days != 1:
            print_msg = ("\tExpt date %s does not match dissection date %s" % (expt_date, dis_time.date()))
            return print_msg
        dis_time_1 = "{d.year}-{d.month}-{d.day} {d.hour}:{d.minute:02d}".format(d=dis_time)
    else:
        dis_time_1 = "{d.hour}:{d.minute:02d}".format(d=dis_time)

    dis_time_2 = dh.info().get('time_of_dissection', '')  # recorded by rig operator
    if dis_time_2 == '':
        if set_data is True:
            dh.setInfo(time_of_dissection=dis_time_1)
            return None
        else:
            print_msg = ("\tSet dissection time: %s" % dis_time_1)
            return print_msg
    else:
        if dis_time_2 != dis_time_1:
            print_msg = ("\tDissection time disagreement:  %r != %r" % (dis_time_2, dis_time_1))
            return print_msg
    return None

def project_check(dh, sub_id, expt_date, species, set_data=False):
    mouse_prod = datetime.date(2017, 10, 01)
    project = dh.info().get('project', None)
    if project is None:
        if species is None and sub_id is None:
            print_msg = ("\tNo specimen, can't set project code")
            return print_msg
        if species.lower() == 'human':
            if set_data is True:
                dh.setInfo(project='human coarse matrix')
                return None
            else:
                print_msg = ("\tSet Project Code: human coarse matrix")
                return print_msg
        if expt_date < mouse_prod:
            if set_data is True:
                dh.setInfo(project='mouse V1 pre-production')
                return None
            else:
                print_msg = ("\tSet Project Code: mouse V1 pre-production")
                return print_msg
        else:
            if set_data is True:
                dh.setInfo(project='mouse V1 coarse matrix')
                return None
            else:
                print_msg = ("\tSet Project Code: mouse V1 coarse matrix")
                return print_msg
    return None

def region_check(dh, species, set_data=False):
    region = dh.info().get('target_region', None)
    if species is None:
        print_msg = ("\tCan't set target region")
        return print_msg
    if species != 'human':
        if region is None:
            if set_data is True:
                dh.setInfo(target_region='V1')
                return None
            else:
                print_msg = ("\tSet target region: V1")
                return print_msg
        if region != 'V1':
            print_msg = ("\tTarget region mismatch: %s != V1" % region)
            return print_msg
    return None

def internal_check(dh, set_data=False):
    internal = dh.info().get('internal', None)
    if internal is not None and internal not in INTERNAL_RECIPES:
        print_msg = ("\tInternal mismatch: %s not in recipe list" % internal)
        return print_msg
    if internal in (None, ''):
        if set_data is True:
            dh.setInfo(internal='Standard K-Gluc')
            return None
        else:
            print_msg = ("\tSet Internal: Standard K-Gluc")
            return print_msg
    return None

def dye_check(dh, species, genotype, set_data=False):
    internal_dye = dh.info().get('internal_dye', None)
    if internal_dye in(None, '') and species is None:
        print_msg = ("\tCan't set internal dye")
        return print_msg
    if internal_dye in(None, '') and species.lower() == 'human':
        if set_data is True:
            dh.setInfo(internal_dye='AF488')
            return None
        else:
            print_msg = ("\tSet internal dye: AF488")
            return print_msg
    if internal_dye in(None, '') and genotype is not None:
        if len(genotype.split(';')) < 3:
            if set_data is True:
                dh.setInfo(internal_dye='AF488')
                return None
            else:
                print_msg = ("\tSet internal dye: AF488, %s looks likes single transgenic" % genotype)
        elif len(genotype.split(';')) >= 3:
            if set_data is True:
                dh.setInfo(internal_dye='Cascade Blue')
                return None
            else:
                print_msg = ("\tSet internal dye: Cascade Blue, %s looks liked quad" % genotype)
        else:
            print_msg = ("\tCan't parse genotype %s, set internal dye manually" % genotype)
        return print_msg
    return None

def rig_check(dh, set_data=False):
    rig = dh.info().get('rig_name', None)
    if rig in (None, ''):
        if set_data is True:
            dh.setInfo(rig_name=config.rig_name)
            return None
        else:
            print_msg = ("\tSet Rig: %s" % config.rig_name)
            return print_msg
    if rig != config.rig_name:
        if set_data is True:
            print_msg = ("\t Rig name mismatch, overrode and set to %s from %s" % (config.rig_name, rig))
            dh.setInfo(rig_name=config.rig_name)
            return print_msg
        else:
            print_msg = ("\tRig mismatch: %s != %s" % (rig, config.rig_name))
            return print_msg
    return None


root = sys.argv[1]
parser = argparse.ArgumentParser()
parser.add_argument('--set-data', action='store_true', default=False, dest='set-data')
args = vars(parser.parse_args(sys.argv[2:]))
set_data = args['set-data']

# find all subject folders that contain at least one site folder
sites = glob.glob(os.path.join(root, '*', 'slice_*', 'site_*'))
checked_days = set()
checked_slices = set()

for path in sites:
    site_dh = getDirHandle(path)
    slice_dh = site_dh.parent()
    day_dh = slice_dh.parent()
    sub_id = day_dh.info().get('animal_ID', None)
    expt_date = datetime.datetime.fromtimestamp(day_dh.info()['__timestamp__']).date()
    species = day_dh.info().get('species', None)
    if species != 'human':
        genotype = species
    if sub_id is not None and species is None:
        try:
            species = day_dh.info().get('LIMS_specimen_info')['organism']
            genotype = day_dh.info().get('LIMS_specimen_info')['genotype']
        except TypeError:
            try:
                species = day_dh.info().get('LIMS_donor_info')['organism']
            except TypeError:
                guess = sub_id[0]
                if guess == 'H':
                    species = 'human'
                else:
                    species = 'mouse'

    header_msg = ("Experiment Date: %s\nAnimal ID: %s" % (expt_date, sub_id))
    print_msg = []
    if day_dh not in checked_days:
        print_msg.append(rig_check(day_dh, set_data=set_data))
        print_msg.append(dissection_check(day_dh, sub_id, expt_date, diss_columns, set_data=set_data))
        print_msg.append(solution_check(day_dh, expt_date, osm_columns, set_data=set_data))
        print_msg.append(region_check(day_dh, species, set_data=set_data))
        print_msg.append(internal_check(day_dh, set_data=set_data))
        print_msg.append(dye_check(day_dh, species, genotype, set_data=set_data))
        checked_days.add(day_dh)

    if slice_dh not in checked_slices:
        print_msg.append(project_check(slice_dh, sub_id, expt_date, species, set_data=set_data))
        checked_slices.add(slice_dh)

    if all([t is None for t in print_msg]):
        continue
    else:
        print header_msg
        for msg in print_msg:
            if msg is not None:
                print msg
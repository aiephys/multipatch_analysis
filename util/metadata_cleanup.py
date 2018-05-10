import os, sys, glob, datetime, csv, re
from acq4.util.DataManager import getDirHandle
from multipatch_analysis.experiment import Experiment

slice_csv_file = 'Slice_Preparation_Record.csv'
osm_csv_file = 'MultiPatch Throughput.xlsx - Solution Making.csv'

# read all dissection times
recs = list(csv.reader(open(slice_csv_file).readlines()))
id_col = recs[0].index('LabTracks ID')
date_col = recs[0].index('Date Sliced')
sac_col = recs[0].index('Time Mouse Sacrificed')
time_col = recs[0].index('heated ACSF1 incub. Comp.')
diss_times = {rec[id_col]: rec for rec in recs[1:]}

# read solution making
osm_recs = list(csv.reader(open(osm_csv_file).readlines()))
osm_date_col = osm_recs[0].index('Date')
recipe_col = osm_recs[0].index('Recipe')
osm_col = osm_recs[0].index('osmolarity')
osm_dates = {osm_rec[osm_date_col]: osm_rec for osm_rec in osm_recs[1:]}

root = sys.argv[1]

# find all subject folders that contain at least one site folder
sites = glob.glob(os.path.join(root, '*', 'slice_*', 'site_*'))
subs = sorted(list(set([os.path.abspath(os.path.join(s, '..', '..')) for s in sites])))
print_msg = ''
for path in subs:
    dh = getDirHandle(path)
    sub_id = dh.info().get('animal_ID', None)
    expt_date = datetime.datetime.fromtimestamp(dh.info()['__timestamp__']).date()
    print("Experiment Date: %s, Animal ID: %s \r" % (expt_date, sub_id))
    if sub_id is None:
        print("\tNo animal_ID for %s\r" % path)
        #continue
    dis_rec = diss_times.get(sub_id, None)
    if dis_rec is None:
        print("\tNo tissue processing record\r")
        #continue
    
    date = dis_rec[date_col]
    sac = dis_rec[sac_col]
    time = dis_rec[time_col]
    sac_time = datetime.datetime.strptime(date + ' ' + sac, '%m/%d/%Y %I:%M %p')
    if time[-2:] in ('AM', 'PM'):
        dis_time = datetime.datetime.strptime(date + ' ' + time, '%m/%d/%Y %I:%M %p')
    else:
        dis_time = datetime.datetime.strptime(date + ' ' + time, '%m/%d/%Y %I:%M')
        if dis_time < sac_time:
            dis_time = dis_time + datetime.timedelta(seconds=12*3600)

    if (dis_time - sac_time).total_seconds() > 3600:
        print("\tSac. / dissection times for are too far apart\r")
        #continue
    if dis_time <= sac_time:
        print("\tDissection before sac\r")
        #continue

    if expt_date != dis_time.date():
        if (expt_date - dis_time.date()).days != 1:
            print("\tExpt date %s does not match dissection date %s\r" % (expt_date, dis_time.date()))
            #continue
        dis_time_1 = "{d.year}-{d.month}-{d.day} {d.hour}:{d.minute:02d}".format(d=dis_time)
    else:
        dis_time_1 = "{d.hour}:{d.minute:02d}".format(d=dis_time)

    dis_time_2 = dh.info().get('time_of_dissection', '')  # recorded by rig operator
    if dis_time_2 == '':
        print("\tSet dissection time: %s\r" % dis_time_1)
        #continue
    else:
        if dis_time_2 != dis_time_1:
            print("\tDissection time disagreement:  %r != %r\r" % (dis_time_2, dis_time_1))
            #continue
        else:
            print("\tDissection Time: Match!\r")
            #continue
    # repeat for solution checking
    solution = dh.info().get('solution', None)
    osm = dh.info().get('solution_osm', None)
    try:
        osm_entry = osm_dates[date]
        recipe = osm_entry[recipe_col]
        sheet_osm = osm_entry[osm_col]
        if solution is None:
            print_msg += ("\tSet aCSF: %s\r" % recipe)
        else:
            calcium = solution.split('m')[0]
            if calcium in recipe:
                print("\taCSF: Match!\r")
            else:
                print("\taCSF disagreement\r")
        if osm is None:
            print("\tSet osmolarity: %s\r" % osmolarity)
        elif sheet_osm == osm:
            print("\tOsmolarity: Match!\r")
        else:
            print("\tOsmolarity disagreement\r")
    except KeyError:
        print("\tDate not in Multipatch_throughput sheet, check notes\r")


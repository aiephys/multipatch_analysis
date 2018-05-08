import os, sys, glob, datetime, csv, re
from acq4.util.DataManager import getDirHandle
from multipatch_analysis.experiment import Experiment

csv_file = 'Slice_Preparation_Record.csv'

# read all dissection times
recs = list(csv.reader(open(csv_file).readlines()))
id_col = recs[0].index('LabTracks ID')
date_col = recs[0].index('Date Sliced')
sac_col = recs[0].index('Time Mouse Sacrificed')
time_col = recs[0].index('heated ACSF1 incub. Comp.')
diss_times = {rec[id_col]: rec for rec in recs[1:]}


root = sys.argv[1]

# find all subject folders that contain at least one site folder
sites = glob.glob(os.path.join(root, '*', 'slice_*', 'site_*'))
subs = sorted(list(set([os.path.abspath(os.path.join(s, '..', '..')) for s in sites])))

for path in subs:
    dh = getDirHandle(path)
    sub_id = dh.info().get('animal_ID', None)
    if sub_id is None:
        print("No animal_ID for %s" % path)
        continue
    dis_rec = diss_times.get(sub_id, None)
    if dis_rec is None:
        print("No tissue processing record for %s" % sub_id)
        continue
    
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
        print("Sac. / dissection times for %s are too far apart" % sub_id)
        continue
    if dis_time <= sac_time:
        print("Dissection before sac: %s" % sub_id)
        continue
    
    expt_date = datetime.datetime.fromtimestamp(dh.info()['__timestamp__']).date()
    if expt_date != dis_time.date():
        if (expt_date - dis_time.date()).days != 1:
            print("Expt date %s does not match dissection date %s for %s" % (expt_date, dis_time.date(), sub_id))
            continue
        dis_time_1 = "{d.year}-{d.month}-{d.day} {d.hour}:{d.minute:02d}".format(d=dis_time)
    else:
        dis_time_1 = "{d.hour}:{d.minute:02d}".format(d=dis_time)

    dis_time_2 = dh.info().get('time_of_dissection', '')  # recorded by rig operator
    if dis_time_2 == '':
        print("Set: %s =>  %s" % (sub_id, dis_time_1))
        continue
    else:
        if dis_time_2 != dis_time_1:
            print("Dissection time disagreement for %s:  %r != %r" % (sub_id, dis_time_2, dis_time_1))
            continue
        else:
            print("Match! %s" % sub_id)
            continue


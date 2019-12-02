import glob, os, argparse, sys, csv
import pandas as pd
from acq4.util.DataManager import getHandle
from aisynphys import config
from aisynphys.util import timestamp_to_datetime
from datetime import datetime, timedelta
from aisynphys.data.pipette_metadata import PipetteMetadata
from collections import OrderedDict

all_paths = glob.glob(os.path.join(config.synphys_data, '*.***'))
nucleus = {'+': 'nucleus_present', '-': 'nucleus_absent', '': None}
daily_report_path = '//allen/programs/celltypes/workgroups/synphys/generated_reports'

def generate_daily_report(day):
    errors = []
    file_name = '%s_mps_Transcriptomics_report.xlsx' % datetime.strftime(day, "%y%m%d")
    file_path = daily_report_path + '/' + file_name
    project_code = '102-01-010-10'
    columns = [
        'Patch Tube Name',
        'Blank Fill Date',
        'Patch Date',
        'Library Prep Day1 Date',
        'Species',
        'Specimen ID',
        'Cell Line',
        'ROI Major',
        'ROI Minor',
        'Comments',
        'Project Code',
        ]

    end = day.date() + timedelta(hours=3)
    start = end - timedelta(hours=24)

    expt_paths = get_expts_in_range(all_paths, start, end)
    if len(expt_paths) == 0:
        print('No experiments found within date range %s - %s' % (datetime.strftime(start,"%m/%d/%Y"), datetime.strftime(end,"%m/%d/%Y")))
        return
    site_paths = [glob.glob(os.path.join(path, 'slice_*', 'site_*')) for path in expt_paths]
    site_paths = [sp for paths in site_paths for sp in paths] #flatten site paths if nested list
    
    row_data = []
    for site in site_paths:
        site_source = open(os.path.join(site, 'sync_source')).read()
        errors.append(site_source)
        site_dh = getHandle(site)
        site_info = site_dh.info()
        slice_info = site_dh.parent().info()
        day_info = site_dh.parent().parent().info()
        pip_meta = PipetteMetadata(site)
        headstages = site_info.get('headstages')
        
        if headstages is None:
            errors.append('\tNo recorded headstages')
            continue
        
        site_errors = {}
        for hs, info in headstages.items():
            if info['Tube ID'] == '':
                continue
            errors.append('\t %s' % hs)
            row = {k: None for k in columns}
            blank_fill_date = slice_info.get('blank_fill_date')
            try:
                datetime.strptime(blank_fill_date, "%m/%d/%Y")
            except ValueError:
                errors.append('\t\tblank fill date has improper format')
                blank_fill_date = ''
            row['Blank Fill Date'] = blank_fill_date
            patch_date = timestamp_to_datetime(day_info.get('__timestamp__'))
            row['Patch Date'] = datetime.strftime(patch_date, "%m/%d/%Y") if isinstance(patch_date, datetime) else None
            species = day_info.get('LIMS_donor_info',{}).get('organism')
            row['Species'] = species.capitalize() if isinstance(species, str) else None
            row['Specimen ID'] = day_info.get('animal_ID')
            row['Cell Line'] = day_info.get('LIMS_donor_info', {}).get('genotype')
            row['ROI Major'] = format_roi_major(day_info.get('target_region'))

            pip = pip_meta.pipettes[hs[-1]]
            tube_id = int(info['Tube ID'].split('_')[2])
            row['tube_id'] = tube_id
            row['Patch Tube Name'] = info['Tube ID']
            row['Comments'] = nucleus[info.get('Nucleus', '')]
            row['ROI Minor'] = format_roi_minor(pip['target_layer'])
            row['Project Code'] = project_code

            for k, v in row.items():
                if v is None and k != 'Library Prep Day1 Date':
                    if k == 'Cell Line' and row['Species'] == 'Human':
                        continue
                    errors.append('\t\t%s has no data' % k)
            row_data.append(row)

    print('\n'.join(errors))
    data_df = pd.DataFrame(row_data)
    data_df.set_index('tube_id', inplace=True)
    data_df.sort_index(inplace=True)
    data_df.to_excel(file_path, index=False)    

def format_roi_major(roi):
    if roi == 'V1':
        return 'VISp'
    else:
        return roi

def format_roi_minor(roi):
    if roi == '2/3':
        return 'L2-3'
    elif roi is not '':
        return 'L%s' % roi

def get_expts_in_range(expts, start, end):
    in_range = []
    for expt in expts:
        expt_date = timestamp_to_datetime(float(os.path.basename(expt))).date()
        if expt_date >= start and expt_date <= end:
            in_range.append(expt)

    return in_range

def valid_date_range(date_range):
    start_date = valid_date(date_range[0])
    end_date = valid_date(date_range[1])
    return [start_date, end_date]


def valid_date(d):
    try:
        return datetime.strptime(d, "%Y%m%d")
    except ValueError:
        msg = "Date '{0}' does not match format YYYYMMDD"
        raise argparse.ArgumentTypeError(msg)

if __name__ == '__main__':
    import pyqtgraph as pg
    pg.dbg()

    parser = argparse.ArgumentParser()
    parser.add_argument('--daily', type=valid_date, nargs='?', const=datetime.now(), help="Default date is today, optionally set the date with the format YYYYMMDD")
    parser.add_argument('--monthly', type=valid_date_range, nargs=2, help="Provide date range to create monthly report as YYYYMMDD YYYYMMDD (start, end)")

    args = parser.parse_args(sys.argv[1:])

    if args.daily is not None:
        day = args.daily
        generate_daily_report(day)

    if args.monthly is not None:
        start_date = args.monthly[0]
        end_date = args.monthly[1]
        generate_monthly_report(start_date, end_date)

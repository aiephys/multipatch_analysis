import glob, os, argparse, sys, csv, re
import pandas as pd
import pyqtgraph as pg
from acq4.util.DataManager import getHandle
from aisynphys import config, lims
from aisynphys.util import timestamp_to_datetime
from datetime import datetime, timedelta
from aisynphys.data.pipette_metadata import PipetteMetadata
from aisynphys.data.slice import Slice
from collections import OrderedDict

all_paths = glob.glob(os.path.join(config.synphys_data, '*.***'))
nucleus = {'+': 'nucleus_present', '-': 'nucleus_absent', '': None}
organism = {'Mus musculus': 'Mouse', 'Homo Sapiens': 'Human'}

def generate_daily_report(day):
    """ Generate a daily PatchSeq report for Kim's team. PatchSeq metadata is collected from the acq4 directories
    for every experiment. Only metadata associated with a Patched Cell Container are processed.
    """

    if day == datetime.today().date():
        day = day - timedelta(hours=24)

    file_name = '%s_mps_Transcriptomics_report.xlsx' % datetime.strftime(day, "%y%m%d")
    file_path = config.patchseq_report_path + '/' + file_name
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

    # collect experiments for the specified day
    expt_paths = get_expts_in_range(all_paths, day, day)
    site_paths = [glob.glob(os.path.join(path, 'slice_*', 'site_*')) for path in expt_paths]
    site_paths = [sp for paths in site_paths for sp in paths] #flatten site paths if nested list
    
    row_data = []
    # look through each site directory
    for site in site_paths:
        if os.path.isdir(site) is False:
            continue
        errors = []
        site_source = open(os.path.join(site, 'sync_source')).read()
        errors.append(site_source)
        site_dh = getHandle(site)
        site_info = site_dh.info()
        slice_info = site_dh.parent().info()
        day_info = site_dh.parent().parent().info()
        pip_meta = PipetteMetadata(site)
        headstages = site_info.get('headstages')
        
        # check to make sure there are recorded headstages and patchseq tubes, else move to next site
        if headstages is None:
            print('%s\tNo recorded headstages' % site_source)
            continue
        tubes = [hs['Tube ID'] for hs in headstages.values()] 
        no_tubes = all([t == '' for t in tubes]) 
        if no_tubes:
            continue

        patch_date_dt = timestamp_to_datetime(day_info.get('__timestamp__'))
        patch_date = datetime.strftime(patch_date_dt, "%m/%d/%Y") if isinstance(patch_date_dt, datetime) else None 
        specimen_id = day_info.get('animal_ID')
        species = lims.specimen_species(slice_info.get('specimen_ID'))
        species = organism.get(species) 
        if species == 'Mouse':
            genotype = day_info.get('LIMS_donor_info', {}).get('genotype')
        else:
            genotype = None
        roi_major = format_roi_major(day_info.get('target_region'))

        blank_fill_date = slice_info.get('blank_fill_date', '')
        try:
            datetime.strptime(blank_fill_date, "%m/%d/%Y")
        except ValueError:
            errors.append('\tblank fill date has improper format')
            blank_fill_date = None
        
        # for headstages that have patchseq tubes log metadata
        for hs, info in headstages.items():
            tube_name, tube_id, msg = parse_tube(info, patch_date_dt)
            if tube_name == None:
                if msg is not None:
                    print('\t\t%s '%hs + msg)
                continue
            row = OrderedDict([k, None] for k in columns)

            pip = pip_meta.pipettes[hs[-1]]
            nucleus_state = nucleus[info.get('Nucleus', '')]
            roi_minor = format_roi_minor(pip['target_layer'])

            row.update({
                'Blank Fill Date': blank_fill_date,
                'Patch Date': patch_date,
                'Specimen ID': specimen_id,
                'Species': species,
                'Cell Line': genotype,
                'Patch Tube Name': tube_name,
                'tube_id': tube_id,
                'Comments': nucleus_state,
                'ROI Major': roi_major,
                'ROI Minor': roi_minor,
                'Project Code': project_code,
            })

            # check that all requried columns are filled in
            for k, v in row.items():
                if v is None and k != 'Library Prep Day1 Date':
                    if k == 'Cell Line' and row['Species'] == 'Human':
                        continue
                    row[k] = 'CHECK DATA'
            row_data.append(row)
        if len(errors) > 1:    
            print('\n'.join(errors))
    
    # convert report to a dataframe and export to excel
    report_df = to_df(row_data, report_type='daily')
    if report_df is not None:
        report_df.to_excel(file_path, index=False)

def generate_monthly_report(start_date, end_date):
    """ Generate a monthly PatchSeq report for Shiny. PatchSeq metadata is collected from the acq4 directories
    for every experiment. Only metadata associated with a Patched Cell Container are processed.
    """

    file_name = '%s_%s_mps_metadata_report.xlsx' % (datetime.strftime(start_date, "%y%m%d"), datetime.strftime(end_date, "%y%m%d"))
    file_path = config.patchseq_report_path + '/' + file_name

    required_cols = {
        'tubeID': 'A',
        'patch.date': 'B',
        'rigOperator': 'C',
        'rigNumber': 'D',
        'Fill.Date': 'E',
        'internalFillDate': 'F',
        'creCell': 'H',
        'manualRoi': 'J',
        'postPatch': 'S',
        'endPipetteR': 'T',
    }

    not_required_cols = {
        'pilotName': 'G',
        'autoRoi': 'I',
        'cell_depth': 'K',
        'sliceHealth': 'L',
        'timeWholeCellStart': 'M',
        'timeExtractionStart': 'N',
        'pressureApplied': 'O',
        'timeExtractionEnd': 'P',
        'retractionPressureApplied': 'Q',
        'timeRetractionEnd': 'R',
    }

    # not all columns are required but they must be in a specified order
    columns = required_cols.copy()
    columns.update(not_required_cols)
    columns = [k for k,v in sorted(columns.items(), key=lambda item: item[1])]

    # collect experiments for the date range provided
    expt_paths = get_expts_in_range(all_paths, start_date, end_date)
    site_paths = [glob.glob(os.path.join(path, 'slice_*', 'site_*')) for path in expt_paths]
    site_paths = [sp for paths in site_paths for sp in paths] #flatten site paths if nested list
    
    row_data = []
    # look through each site directory for patchseq data
    for site in site_paths:
        errors = []
        site_source = open(os.path.join(site, 'sync_source')).read()
        errors.append(site_source)
        site_dh = getHandle(site)
        site_info = site_dh.info()
        slice_info = site_dh.parent().info()
        day_dh = site_dh.parent().parent()
        day_info = day_dh.info()
        pip_meta = PipetteMetadata(site)
        headstages = site_info.get('headstages')
        
        # if no headstages were recorded or tubes collected, move along
        if headstages is None:
            print('%s\tNo recorded headstages' % site_source)
            continue
        tubes = [hs['Tube ID'] for hs in headstages.values()] 
        no_tubes = all([t == '' for t in tubes]) 
        if no_tubes:
            continue

        index_file = pg.configfile.readConfigFile(os.path.join(day_dh.path, '.index'))
        rig_name = index_file['.'].get('rig_name')
        patch_date_dt = timestamp_to_datetime(day_info.get('__timestamp__'))
        patch_date = datetime.strftime(patch_date_dt, "%m/%d/%Y") if isinstance(patch_date_dt, datetime) else None
        operator = day_info.get('rig_operator', '')
        roi = format_roi_major(day_info.get('target_region'))
        slic = Slice(site_dh.parent().name())
        genotype = slic.genotype
        if genotype is None and slic.species == 'Mouse':
            errors.append('\tno genotype for %s, this may affect the creCell column' % slic.lims_specimen_name)

        blank_fill_date = slice_info.get('blank_fill_date', '')
        try:
            datetime.strptime(blank_fill_date, "%m/%d/%Y")
        except ValueError:
            errors.append('\tblank fill date has improper format')
            blank_fill_date = None
        
        internal_fill_date = slice_info.get('internal_fill_date', '')
        try:
            datetime.strptime(internal_fill_date, "%m/%d/%Y")
        except ValueError:
            errors.append('\tinternal fill date has improper format')
            internal_fill_date = None
        
        for hs, info in headstages.items():
            tube_name, tube_id, msg = parse_tube(info, patch_date_dt)
            if tube_name == None:
                if msg is not None:
                    print('\t\t%s '%hs + msg)
                continue
            row = OrderedDict([k, None] for k in columns)
            
            human_culture = True if tube_name[1] == 'T' else False
            if human_culture is True and genotype is None:
                errors.append('\tno genotype for %s, this may affect the creCell column' % slic.lims_specimen_name)
            color = info.get('Reporter')
            reporter = None
            if color == '-':
                reporter = color
            elif color in ['red', 'green', 'yellow'] and genotype is not None:
                reporter = genotype.color_to_reporter(color)
            elif color == 'NA':
                reporter = ''
            
            pip = pip_meta.pipettes[hs[-1]]
            layer = pip['target_layer']
            manual_roi = roi + layer if (roi not in [None, ''] and layer not in [None, '']) else None
            nucleus_state = nucleus[info.get('Nucleus', '')]
            end_seal = info['End Seal']
            end_seal = 1000 if end_seal else 0 # in MOhms

            row.update({
                'internalFillDate': internal_fill_date,
                'Fill.Date': blank_fill_date,
                'tubeID': tube_name,
                'tube_id': tube_id,
                'patch.date': patch_date,
                'rigOperator': operator,
                'rigNumber': rig_name,
                'creCell': reporter,
                'manualRoi': manual_roi,
                'postPatch': nucleus_state,
                'endPipetteR': end_seal,
            })

            # check that there is metadata for all required columns
            for k in required_cols.keys():
                v = row[k]
                if v is None:
                    row[k] = 'CHECK DATA'
                    errors.append('\t\t%s %s has no data' % (hs, k))
            row_data.append(row)
        if len(errors) > 1:
            print('\n'.join(errors))

    report_df = to_df(row_data, report_type='monthly')
    
    # cross-check with daily reports to make sure all tubes are accounted for
    tube_cross_check(report_df['tubeID'], start_date, end_date)

    report_df.to_excel(file_path, index=False)  

def tube_cross_check(monthly_tubes, start_date, end_date):
    daily_report_folder = '//allen/programs/celltypes/workgroups/synphys/MPS_transcriptomics_report'
    delta = end_date - start_date
    days = [start_date + timedelta(days=i) for i in range(delta.days + 1)] 
    daily_tubes = None
    for day in days:
        daily_report_file = os.path.join(daily_report_folder, '%s_mps_Transcriptomics_report.xlsx' % datetime.strftime(day, "%y%m%d"))
        try:
            daily_df = pd.read_excel(daily_report_file, header=0, index_col=None)
        except:
            print('No report for %s' % datetime.strftime(day, "%y%m%d"))
            continue
        if daily_tubes is None:
            daily_tubes = daily_df['Patch Tube Name']
        else:
            daily_tubes = daily_tubes.append(daily_df['Patch Tube Name'], ignore_index=True)
    
    daily_tubes = daily_tubes.replace('\n', '', regex=True)
    monthly_tubes = monthly_tubes.replace('\n', '', regex=True)
    extra = monthly_tubes[~monthly_tubes.isin(daily_tubes)]
    missing = daily_tubes[~daily_tubes.isin(monthly_tubes)]

    print('### These tubes appear to be missing from the Monthly Report:\n%s' % missing.to_string(index=False))
    print('### These tubes appear in the Monthly Report but not in a Daily Report:\n%s' % extra.to_string(index=False))

def parse_tube(info, date):
    if info['Tube ID'] == '':
        return None, None, None
    tube_name = info['Tube ID']
    name_check = re.match(r'P(M|T)S4_(?P<date>\d{6})_(?P<tube_id>\d{3})_A01', tube_name)
    if name_check is None:
        msg = 'Tube %s ID does not match the proper format' % tube_name
        return None, None, msg
    date_check = name_check.group('date') == datetime.strftime(date, "%y%m%d")
    if date_check is False:
        msg = 'Date %s in Tube ID does not match the date of the experiment' % name_check.group('date')
        return None, None, msg

    tube_id = name_check.group('tube_id')
    return tube_name, tube_id, None

def to_df(report_data, report_type):
    if len(report_data) == 0:
        print('No patchseq tubes to report')
        if report_type == 'monthly':
            return pd.DataFrame(columns=['tubeID'])
        else:
            return None
    data_df = pd.DataFrame(report_data)
    if report_type == 'daily':
        data_df.sort_values('tube_id', inplace=True)
    elif report_type == 'monthly':
        data_df.sort_values(['patch.date', 'tube_id'], inplace=True)
    data_df.drop('tube_id', axis='columns', inplace=True)  

    return data_df

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

    if len(in_range) == 0:
        print('No experiments found within date range %s - %s' % (datetime.strftime(start,"%m/%d/%Y"), datetime.strftime(end,"%m/%d/%Y")))
        exit()
    return in_range

def valid_date(d):
    try:
        return datetime.strptime(d, "%Y%m%d").date()
    except ValueError:
        msg = "Date '{0}' does not match format YYYYMMDD"
        raise argparse.ArgumentTypeError(msg)

if __name__ == '__main__':
    pg.dbg()

    parser = argparse.ArgumentParser()
    parser.add_argument('--daily', type=valid_date, nargs='?', const=datetime.now().date(), help="Default date is today, optionally set the date with the format YYYYMMDD")
    parser.add_argument('--monthly', type=valid_date, nargs=2, help="Provide date range to create monthly report as YYYYMMDD YYYYMMDD (start, end)")

    args = parser.parse_args(sys.argv[1:])

    if args.daily is not None:
        day = args.daily
        generate_daily_report(day)

    if args.monthly is not None:
        start_date = args.monthly[0]
        end_date = args.monthly[1]
        generate_monthly_report(start_date, end_date)

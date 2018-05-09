from __future__ import print_function
import os, sys, datetime, re, glob, traceback, time, atexit, threading
try:
    import queue
except ImportError:
    import Queue as queue
from pprint import pprint
from collections import OrderedDict
import numpy as np
from multipatch_analysis import config, lims
from multipatch_analysis.experiment import Experiment
from multipatch_analysis.database import database
from multipatch_analysis.genotypes import Genotype
from multipatch_analysis.ui.actions import ExperimentActions

from acq4.util.DataManager import getDirHandle
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore


fail_color = (255, 200, 200)
pass_color = (200, 255, 200)


class Dashboard(QtGui.QWidget):
    def __init__(self, limit=0, no_thread=False, filter_defaults=None):
        QtGui.QWidget.__init__(self)

        # fields displayed in ui
        self.visible_fields = [
            ('timestamp', float), 
            ('rig', 'U100'), 
            ('path', 'U100'), 
            ('project', 'U100'),
            ('description', 'U100'), 
            ('primary', 'U100'), 
            ('archive', 'U100'), 
            ('backup', 'U100'), 
            ('NAS', 'U100'), 
            ('data', 'U100'),
            ('submitted', 'U100'),
            ('connections', 'U100'),
            ('site.mosaic', 'U100'), 
            ('DB', 'U100'), 
            ('LIMS', 'U100'), 
            ('20x', 'U100'), 
            ('cell map', 'U100'), 
            ('63x', 'U100'), 
            ('morphology', 'U100')
        ]

        # data tracked but not displayed
        self.hidden_fields = [
            ('experiment', object),
            ('item', object),
            ('error', object),
        ]

        # maps field name : index (column number)
        self.field_indices = {self.visible_fields[i][0]:i for i in range(len(self.visible_fields))}

        self.records = GrowingArray(dtype=self.visible_fields + self.hidden_fields)
        self.records_by_expt = {}  # maps expt:index

        self.selected = None

        # set up UI
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.splitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        self.layout.addWidget(self.splitter, 0, 0)

        self.filter = pg.DataFilterWidget()
        self.filter_fields = OrderedDict([(field[0], {'mode': 'enum', 'values': {}}) for field in self.visible_fields])
        self.filter_fields['timestamp'] = {'mode': 'range'}
        self.filter_fields.pop('path')
        if filter_defaults is not None:
            for k, defs in filter_defaults.items():
                self.filter_fields[k]['values'].update(defs)

        self.filter.setFields(self.filter_fields)
        self.splitter.addWidget(self.filter)
        self.filter.sigFilterChanged.connect(self.filter_changed)
        self.filter.addFilter('rig')
        self.filter.addFilter('project')

        self.expt_tree = pg.TreeWidget()
        self.expt_tree.setSortingEnabled(True)
        self.expt_tree.setColumnCount(len(self.visible_fields))
        self.expt_tree.setHeaderLabels([f[0] for f in self.visible_fields])
        self.splitter.addWidget(self.expt_tree)
        self.expt_tree.itemSelectionChanged.connect(self.tree_selection_changed)

        self.status = QtGui.QStatusBar()
        self.status.setMaximumHeight(25)
        self.layout.addWidget(self.status, 1, 0)

        self.resize(1000, 900)
        self.splitter.setSizes([200, 800])
        colsizes = [150, 50, 230, 150, 200]
        for i in range(len(self.visible_fields)):
            size = colsizes[i] if i < len(colsizes) else 50
            self.expt_tree.setColumnWidth(i, size)

        # Add reload context menu action
        self.setContextMenuPolicy(QtCore.Qt.ActionsContextMenu)
        self.reload_action = QtGui.QAction('Reload', self)
        self.reload_action.triggered.connect(self.reload_clicked)
        self.addAction(self.reload_action)
        
        # Add generic experiment context menu actions
        self.expt_actions = ExperimentActions()
        for act in self.expt_actions.actions.values():
            self.addAction(act)
        
        # Queue of experiments to be checked
        self.expt_queue = queue.PriorityQueue()

        # collect a list of all data sources to search
        search_paths = [config.synphys_data]
        for rig_name, rig_path_sets in config.rig_data_paths.items():
            for path_set in rig_path_sets:
                search_paths.extend(list(path_set.values()))

        # Each poller searches a single data source for new experiments, adding them to the queue
        self.pollers = []
        self.poller_status = {}
        for search_path in search_paths:
            if not os.path.exists(search_path):
                print("Ignoring search path:", search_path)
                continue
            poll_thread = PollThread(self.expt_queue, search_path, limit=limit)
            poll_thread.update.connect(self.poller_update)
            if no_thread:
                poll_thread.poll()  # for local debugging
            else:
                poll_thread.start()
            self.pollers.append(poll_thread)

        # Checkers pull experiments off of the queue and check their status
        if no_thread:
            self.checker = ExptCheckerThread(self.expt_queue)
            self.checker.update.connect(self.checker_update)
            self.checker.run(block=False)            
        else:
            self.checkers = []
            for i in range(3):
                self.checkers.append(ExptCheckerThread(self.expt_queue))
                self.checkers[-1].update.connect(self.checker_update)
                self.checkers[-1].start()

        # shut down threads nicely on exit
        atexit.register(self.quit)

        self.status_timer = QtCore.QTimer()
        self.status_timer.timeout.connect(self.update_status)
        self.status_timer.start(1000)

    def tree_selection_changed(self):
        sel = self.expt_tree.selectedItems()[0]
        rec = self.records[sel.index]
        self.selected = rec
        expt = rec['experiment']
        self.expt_actions.experiment = expt
        print("===================", expt.timestamp)
        print(" description:", rec['description'])
        print("   spec name:", expt.specimen_name)
        print("    genotype:", expt.lims_record.get('genotype', None))
        print("    NAS path:", expt.nas_path)
        print("primary path:", expt.primary_path)
        print("archive path:", expt.archive_path)
        print(" backup path:", expt.backup_path)
        err = rec['error']
        if err is not None:
            print("Error checking experiment:")
            traceback.print_exception(*err)

    def poller_update(self, path, status):
        """Received an update on the status of a poller
        """
        self.poller_status[path] = status
        paths = sorted(self.poller_status.keys())
        msg = '    '.join(['%s: %s' % (p,self.poller_status[p]) for p in paths])
        # self.long_status.showMessage(msg)

    def update_status(self):
        self.status.showMessage("%d sites in queue" % self.expt_queue.qsize())

    def checker_update(self, rec):
        """Received an update from a worker thread describing information about an experiment
        """
        expt = rec['experiment']
        if expt in self.records_by_expt:
            # use old record / item
            index = self.records_by_expt[expt]
            item = self.records[index]['item']
        else:
            # add new record / item
            item = pg.TreeWidgetItem()
            self.expt_tree.addTopLevelItem(item)
            rec['item'] = item
            index = self.records.add_record({})
            item.index = index

        record = self.records[index]

        # update item/record fields
        update_filter = False
        for field, val in rec.items():
            if field in self.field_indices and isinstance(val, tuple):
                # if a tuple was given, interpret it as (text, color)
                val, color = val
            else:
                # otherwise make a guess on a good color
                color = None
                if val is True:
                    color = pass_color
                elif val is False:
                    color = fail_color
                elif val in ('ERROR', 'MISSING'):
                    color = fail_color
            display_val = str(val)

            if field == 'timestamp' and rec.get('error') is not None:
                color = fail_color

            # update this field in the record
            record[field] = val
            
            # update this field in the tree item
            try:
                i = self.field_indices[field]
            except KeyError:
                continue

            item.setText(i, display_val)
            if color is not None:
                item.setBackgroundColor(i, pg.mkColor(color))

            # update filter fields
            filter_field = self.filter_fields.get(field)
            if filter_field is not None and filter_field['mode'] == 'enum' and val not in filter_field['values']:
                filter_field['values'][val] = True
                update_filter = True

        if update_filter:
            self.filter.setFields(self.filter_fields)

    def filter_changed(self):
        mask = self.filter.generateMask(self.records)
        for i,item in enumerate(self.records['item']):
            item.setHidden(not mask[i])

    def quit(self):
        print("Stopping pollers..")
        for t in self.pollers:
            t.stop()
            t.wait()

        print("Stopping checkers..")
        for t in self.checkers:
            t.stop()
        for t in self.checkers:
            t.wait()  # don't wait until all checkers have been requested to stop!

    def reload_clicked(self, *args):
        expt = self.selected['experiment']
        self.expt_queue.put((-expt.timestamp, expt))


class GrowingArray(object):
    def __init__(self, dtype, init_size=1000):
        self.size = 0
        self._data = np.empty(init_size, dtype=dtype)
        self._view = self._data[:0]

    def __len__(self):
        return self.size

    @property
    def shape(self):
        return (self.size,)

    def __getitem__(self, i):
        return self._view[i]

    def add_record(self, rec):
        index = self.size
        self._grow(self.size+1)
        self.update_record(index, rec)
        return index

    def update_record(self, index, rec):
        for k,v in rec.items():
            print(k, v)
            self._data[index][k] = v

    def _grow(self, size):
        if size > len(self._data):
            self._data = np.resize(self._data, len(self._data)*2)
        self._view = self._data[:size]
        self.size = size


class PollThread(QtCore.QThread):
    """Used to check in the background for changes to experiment status.
    """
    update = QtCore.Signal(object, object)  # search_path, status_message
    known_expts = {}
    known_expts_lock = threading.Lock()

    def __init__(self, expt_queue, search_path, limit=0):
        QtCore.QThread.__init__(self)
        self.expt_queue = expt_queue
        self.search_path = search_path
        self.limit = limit
        self._stop = False
        self.waker = threading.Event()
        
    def stop(self):
        self._stop = True
        self.waker.set()

    def run(self):
        while True:
            try:
                # check for new experiments hourly
                self.poll()
                self.waker.wait(3600)
                if self._stop:
                    return
            except Exception:
                sys.excepthook(*sys.exc_info())
            break
                
    def poll(self):
        expts = {}

        # Find all available site paths across all data sources
        count = 0
        path = self.search_path

        self.update.emit(path, "Updating...")
        root_dh = getDirHandle(path)

        # iterate over all expt sites in this path
        for expt_path in glob.iglob(os.path.join(root_dh.name(), '*', 'slice_*', 'site_*')):
            if self._stop:
                return

            expt = ExperimentMetadata(path=expt_path)
            ts = expt.timestamp

            # Couldn't get timestamp; show an error message
            if ts is None:
                print("Error getting timestamp for %s" % expt)
                continue

            with self.known_expts_lock:
                if ts in self.known_expts:
                    # We've already seen this expt elsewhere; skip
                    continue            
                self.known_expts[ts] = expt

            # Add this expt to the queue to be checked
            self.expt_queue.put((-ts, expt))

            count += 1
            if self.limit > 0 and count >= self.limit:
                return
        self.update.emit(path, "Finished")


class ExptCheckerThread(QtCore.QThread):
    update = QtCore.Signal(object)

    def __init__(self, expt_queue):
        QtCore.QThread.__init__(self)
        self.expt_queue = expt_queue
        self._stop = False

    def stop(self):
        self._stop = True
        self.expt_queue.put(('stop', None))

    def run(self, block=True):
        while True:
            ts, expt = self.expt_queue.get(block=block)
            if self._stop or ts == 'stop':
                return
            rec = self.check(expt)
            self.update.emit(rec)

    def check(self, expt):
        try:
            rec = {'experiment': expt}

            rec['timestamp'] = expt.timestamp
            rec['rig'] = expt.rig_name
            rec['path'] = expt.expt_subpath
            rec['project'] = expt.slice_info.get('project', None)
            rec['primary'] = False if expt.primary_path is None else (True if os.path.exists(expt.primary_path) else "-")
            rec['archive'] = False if expt.archive_path is None else (True if os.path.exists(expt.archive_path) else "MISSING")
            rec['backup'] = False if expt.backup_path is None else (True if os.path.exists(expt.backup_path) else "MISSING")
            rec['NAS'] = False if expt.nas_path is None else (True if os.path.exists(expt.nas_path) else "MISSING")

            org = expt.organism
            if org is None:
                description = ("no LIMS spec info", fail_color)
            elif org == 'human':
                description = org
            else:
                try:
                    gtyp = expt.genotype
                    description = ','.join(gtyp.drivers())
                except Exception:
                    description = (org + ' (unknown: %s)'%expt.lims_record.get('genotype', None), fail_color)
            rec['description'] = description
            
            search_paths = [expt.primary_path, expt.archive_path, expt.nas_path]
            submitted = False
            for path in search_paths:
                if path is None:
                    continue
                if os.path.isfile(os.path.join(path, 'pipettes.yml')):
                    submitted = True
                    break

            rec['submitted'] = submitted
            rec['data'] = '-' if expt.nwb_file is None else True
            if rec['submitted']:
                if rec['data'] is True:
                    rec['site.mosaic'] = expt.mosaic_file is not None            
                    rec['DB'] = expt.in_database

                    subs = expt.lims_submissions
                    if subs is None:
                        lims = "ERROR"
                    else:
                        lims = len(subs) == 1
                    rec['LIMS'] = lims
        
        except Exception as exc:
            rec['error'] = sys.exc_info()
        
        return rec


class ExperimentMetadata(Experiment):
    """Handles reading experiment metadata from several possible locations.
    """
    def __init__(self, path=None):
        Experiment.__init__(self, verify=False)
        self._site_path = path

        self.site_dh = getDirHandle(path)
        # self.site_info = self.site_dh.info()
        # self._slice_info = None
        # self._expt_info = None
        # self._specimen_info = None
        self._rig_name = None
        self._primary_path = None
        self._archive_path = None
        self._backup_path = None

    def _get_raw_paths(self):
        expt_subpath = self.expt_subpath

        # find the local primary/archive paths that contain this experiment
        found_paths = False
        rig_name = self.rig_name
        if rig_name is None:
            return
        rig_data_paths = config.rig_data_paths.get(rig_name, [])
        for path_set in rig_data_paths:
            for root in path_set.values():
                test_path = os.path.join(root, expt_subpath)
                if not os.path.isdir(test_path):
                    continue
                dh = getDirHandle(test_path)
                if dh.info()['__timestamp__'] == self.site_info['__timestamp__']:
                    found_paths = True
                    # set self._primary_path, self._archive_path, etc.
                    for k,v in path_set.items():
                        setattr(self, '_'+k+'_path', os.path.join(v, expt_subpath))
                    break
            if found_paths:
                break

    @property
    def expt_subpath(self):
        path = self.site_dh.name()

        # get raw data subpath  (eg: 2018.01.20_000/slice_000/site_000)
        if os.path.abspath(path).startswith(os.path.abspath(config.synphys_data).rstrip(os.path.sep) + os.path.sep):
            # this is a server path; need to back-convert to rig path
            source_path = open(os.path.join(path, 'sync_source')).read()
            expt_subpath = os.path.join(*source_path.split('/')[-3:])
            assert os.path.abspath(path) == os.path.abspath(self.nas_path), "Expected equal paths:\n  %s\n  %s" % (os.path.abspath(path), os.path.abspath(self.nas_path))
        else:
            expt_subpath = os.path.join(*path.split(os.path.sep)[-3:])

        return expt_subpath

    @property
    def primary_path(self):
        if self._primary_path is None:
            self._get_raw_paths()
        return self._primary_path

    @property
    def archive_path(self):
        if self._archive_path is None:
            self._get_raw_paths()
        return self._archive_path

    @property
    def backup_path(self):
        if self._backup_path is None:
            self._get_raw_paths()
        return self._backup_path

    @property
    def slice_dh(self):
        return self.site_dh.parent()

    @property
    def expt_dh(self):
        return self.slice_dh.parent()

    @property
    def nas_path(self):
        expt_dir = '%0.3f' % self.expt_info['__timestamp__']
        subpath = self.site_dh.name(relativeTo=self.expt_dh)
        return os.path.abspath(os.path.join(config.synphys_data, expt_dir, subpath))

    @property
    def organism(self):
        org = self.expt_info.get('organism')
        if org is not None:
            return org
        spec_info = self.lims_record
        if spec_info is None:
            return None
        return spec_info['organism']

    @property
    def rig_name(self):
        if self._rig_name is None:
            name = self.expt_info.get('rig_name')
            if name is not None:
                self._rig_name = name.lower()
            else:
                # infer rig name from paths
                nas_path = self.nas_path
                if nas_path is not None:
                    sync_file = os.path.join(self.nas_path, 'sync_source')
                else:
                    sync_file = None
                if sync_file is not None and os.path.isfile(sync_file):
                    path = open(sync_file, 'rb').read()
                else:
                    path = self.site_dh.name()
                m = re.search(r'(/|\\)(mp\d)(/|\\)', path)
                if m is None:
                    path = os.path.abspath(self.site_dh.name())
                    for rig,paths in config.rig_data_paths.items():
                        data_paths = paths.values()
                        for data_path in paths:
                            if path.startswith(os.path.abspath(data_path)):
                                self._rig_name = rig
                                break
                        if self._rig_name is not None:
                            break
                    
                    return None
                self._rig_name = m.groups()[1].lower()
        return self._rig_name

    @property
    def in_database(self):
        session = database.Session()
        expts = session.query(database.Experiment.id).filter(database.Experiment.acq_timestamp==self.datetime).all()
        session.close()
        return len(expts) == 1

    @property
    def lims_submissions(self):
        if self.lims_record is None:
            return None
        spec_id = self.lims_record['specimen_id']
        if spec_id is None:
            return None
        return lims.expt_submissions(spec_id, self.timestamp)

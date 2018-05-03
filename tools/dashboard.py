import os, sys, datetime, re, glob, traceback, time, queue
from pprint import pprint
from multipatch_analysis import config, lims
from multipatch_analysis.database import database
from multipatch_analysis.genotypes import Genotype
from acq4.util.DataManager import getDirHandle
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore


fail_color = (255, 200, 200)
pass_color = (200, 255, 200)


class Dashboard(QtGui.QWidget):
    def __init__(self, limit=0, no_thread=False):
        QtGui.QWidget.__init__(self)
        
        self.records = {}

        # set up UI
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.expt_tree = pg.TreeWidget()
        self.expt_tree.setSortingEnabled(True)
        self.fields = ['timestamp', 'path', 'rig', 'description', 'primary', 'archive', 'backup', 'NAS', 'pipettes.yml', 'site.mosaic', 'DB', 'LIMS', '20x', 'cell map', '63x', 'morphology']
        self.expt_tree.setColumnCount(len(self.fields))
        self.expt_tree.setHeaderLabels(self.fields)
        self.layout.addWidget(self.expt_tree, 0, 0)
        self.expt_tree.itemSelectionChanged.connect(self.tree_selection_changed)
        
        self.resize(1000, 900)
        
        # Queue of experiments to be checked
        self.expt_queue = queue.PriorityQueue()

        # Poller searches all data sources for new experiments, adding them to the queue
        self.poll_thread = PollThread(self.expt_queue, limit=limit)
        self.poll_thread.update.connect(self.poller_update)
        if no_thread:
            self.poll_thread.poll()  # for local debugging
        else:
            self.poll_thread.start()

        # Checkers pull experiments off of the queue and check their status
        if no_thread:
            self.checker = ExptCheckerThread(self.expt_queue)
            self.checker.update.connect(self.poller_update)
            self.checker.run(block=False)
            
        else:
            self.checker_threads = []
            for i in range(3):
                self.checker_threads.append(ExptCheckerThread(self.expt_queue))
                self.checker_threads[-1].update.connect(self.poller_update)
                self.checker_threads[-1].start()

    def tree_selection_changed(self):
        sel = self.expt_tree.selectedItems()[0]
        pprint(sel.rec)

    def poller_update(self, rec):
        """Received an update from a worker thread describing information about an experiment
        """
        expt = rec['expt']
        if expt in self.records:
            item = self.records[expt]['item']
            self.records[expt].update(rec)
        else:
            item = pg.TreeWidgetItem()
            self.expt_tree.addTopLevelItem(item)
            rec['item'] = item
            self.records[expt] = rec
            item.expt = rec['expt']
            item.rec = rec
            
        for field, val in rec.items():
            try:
                i = self.fields.index(field)
            except ValueError:
                continue
            if isinstance(val, tuple):
                val, color = val
            else:
                color = None
                if val is True:
                    color = pass_color
                elif val in (False, 'ERROR', 'MISSING'):
                    color = fail_color

            item.setText(i, str(val))
            if color is not None:
                item.setBackgroundColor(i, pg.mkColor(color))


class PollThread(QtCore.QThread):
    """Used to check in the background for changes to experiment status.
    """
    update = QtCore.Signal(object)
    
    def __init__(self, expt_queue, limit=0):
        self.expt_queue = expt_queue
        self.limit = limit
        QtCore.QThread.__init__(self)
        
    def run(self):
        while True:
            try:
                # check for new experiments hourly
                self.poll()
                time.sleep(3600)
            except Exception:
                sys.excepthook(*sys.exc_info())
            break
                
    def poll(self):

        self.session = database.Session()
        expts = {}

        print("loading site paths..")
        # collect a list of all data sources to search
        search_paths = [config.synphys_data]
        for rig_name, rig_path_sets in config.rig_data_paths.items():
            for path_set in rig_path_sets:
                search_paths.extend(list(path_set.values()))

        # Find all available site paths across all data sources
        count = 0
        for path in search_paths:
            print("====  Searching %s" % path)
            if not os.path.exists(path):
                print("    skipping; path does not exist.")
                continue
            root_dh = getDirHandle(path)

            # iterate over all expt sites in this path
            for expt_path in glob.iglob(os.path.join(root_dh.name(), '*', 'slice_*', 'site_*')):

                expt = ExperimentMetadata(path=expt_path)
                ts = expt.timestamp

                # Couldn't get timestamp; show an error message
                if ts is None:
                    print("Error getting timestamp for %s" % expt)
                    continue

                # We've already seen this expt elsewhere; skip
                if ts in expts:
                    continue
                
                # Add this expt to the queue
                expts[ts] = expt
                self.expt_queue.put((-ts, expt))

                count += 1
                if self.limit > 0 and count >= self.limit:
                    print("Hit limit; exiting poller")
                    return


class ExptCheckerThread(QtCore.QThread):
    update = QtCore.Signal(object)

    def __init__(self, expt_queue):
        QtCore.QThread.__init__(self)
        self.expt_queue = expt_queue

    def run(self, block=True):
        while True:
            ts, expt = self.expt_queue.get(block=block)
            try:
                rec = self.check(expt)
                self.update.emit(rec)
            except Exception:
                rec = {
                    'expt': expt,
                    'path': expt.site_dh.name(),
                    'timestamp': expt.timestamp or 'ERROR',
                    'rig': 'ERROR',
                }
                self.update.emit(rec)
                traceback.print_exc()
                print(expt)

    def check(self, expt):
        print("   check %s" % expt.site_dh.name())

        org = expt.organism
        if org is None:
            description = ("no LIMS spec info", fail_color)
        elif org == 'human':
            description = org
        else:
            gtyp = expt.genotype
            if gtyp is None:
                description = (org + ' (no genotype)', fail_color)
            else:
                try:
                    description = ','.join(Genotype(gtyp).drivers())
                except Exception:
                    description = (org + ' (unknown: %s)'%gtyp, fail_color)

        subs = expt.lims_submissions
        if subs is None:
            lims = "ERROR"
        else:
            lims = len(subs) == 1

        rec = {
            'expt': expt,
            'path': expt.site_dh.name(), 
            'timestamp': expt.timestamp, 
            'rig': expt.rig_name, 
            'primary': False if expt.primary_path is None else (True if os.path.exists(expt.primary_path) else "MISSING"),
            'archive': False if expt.archive_path is None else (True if os.path.exists(expt.archive_path) else "MISSING"),
            'backup': False if expt.backup_path is None else (True if os.path.exists(expt.backup_path) else "MISSING"),
            'description': description,
            'pipettes.yml': expt.pipette_file is not None,
            'site.mosaic': expt.mosaic_file is not None,
            'DB': expt.in_database,
            'NAS': expt.nas_path is not None,
            'LIMS': lims,
        }
        return rec


class ExperimentMetadata(object):
    """Handles reading experiment metadata from several possible locations.
    """
    def __init__(self, path=None):
        self.site_dh = getDirHandle(path)
        self.site_info = self.site_dh.info()
        self._slice_info = None
        self._expt_info = None
        self._specimen_info = None
        self._rig_name = None
        self._primary_path = None
        self._archive_path = None
        self._backup_path = None

    def _get_raw_paths(self):
        path = self.site_dh.name()

        # get raw data subpath  (eg: 2018.01.20_000/slice_000/site_000)
        if os.path.abspath(path).startswith(os.path.abspath(config.synphys_data) + os.path.sep):
            # this is a server path; need to back-convert to rig path
            source_path = open(os.path.join(path, 'sync_source')).read()
            expt_subpath = os.path.join(*source_path.split('/')[-3:])
            assert os.path.abspath(path) == os.path.abspath(self.nas_path), "Expected equal paths:\n  %s\n  %s" % (os.path.abspath(path), os.path.abspath(self.nas_path))
        else:
            expt_subpath = os.path.join(*path.split(os.path.sep)[-3:])
        
        # find the local primary/archive paths that contain this experiment
        found_paths = False
        rig_data_paths = config.rig_data_paths[self.rig_name]
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
    def slice_info(self):
        if self._slice_info is None:
            self._slice_info = self.slice_dh.info()
        return self._slice_info

    @property
    def expt_info(self):
        if self._expt_info is None:
            self._expt_info = self.expt_dh.info()
        return self._expt_info

    @property
    def specimen_info(self):
        if self._specimen_info is None:
            spec_name = self.slice_info['specimen_ID'].strip()
            try:
                self._specimen_info = lims.specimen_info(spec_name)
            except Exception as exc:
                if 'returned 0 results' in exc.args[0]:
                    pass
                else:
                    raise
        return self._specimen_info

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
        spec_info = self.specimen_info
        if spec_info is None:
            return None
        return spec_info['organism']

    @property
    def genotype(self):
        gtyp = self.expt_info.get('genotype')
        if gtyp is not None:
            return gtyp
        spec_info = self.specimen_info
        if spec_info is None:
            return None
        return spec_info['genotype']

    @property
    def rig_name(self):
        if self._rig_name is None:
            name = self.expt_info.get('rig_name')
            if name is not None:
                self._rig_name = name.lower()
            else:
                # infer rig name from paths
                if 'sync_source' in self.site_dh.ls():
                    path = self.site_dh['sync_source'].read()
                else:
                    path = self.site_dh.name()
                m = re.search(os.path.sep + '(mp\d)', path)
                if m is None:
                    raise Exception("Can't determine rig name for %s" % path)
                self._rig_name = m.groups()[0].lower()
        return self._rig_name

    @property
    def timestamp(self):
        return self.site_info.get('__timestamp__')

    @property
    def datetime(self):
        return datetime.datetime.fromtimestamp(self.timestamp)

    @property
    def pipette_file(self):
        if 'pipettes.yml' in self.site_dh.ls():
            return self.site_dh['pipettes.yml'].name()
        return None

    @property
    def mosaic_file(self):
        if 'site.mosaic' in self.site_dh.ls():
            return self.site_dh['site.mosaic'].name()
        return None

    @property
    def in_database(self):
        session = database.Session()
        expts = session.query(database.Experiment).filter(database.Experiment.acq_timestamp==self.datetime).all()
        return len(expts) == 1

    @property
    def lims_submissions(self):
        if self.specimen_info is None:
            return None
        spec_id = self.specimen_info['specimen_id']
        if spec_id is None:
            return None
        return lims.expt_submissions(spec_id, self.timestamp)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--no-thread', action='store_true', default=False, dest='no_thread',
                    help='Do all polling in main thread (to make debugging easier).')
    parser.add_argument('--limit', type=int, dest='limit', default=0, help="Limit the number of experiments to poll (to make testing easier).")
    args = parser.parse_args(sys.argv[1:])

    app = pg.mkQApp()
    console = pg.dbg()
    db = Dashboard(limit=args.limit, no_thread=args.no_thread)
    db.show()
    
    def sel_change():
        global sel
        sel = db.expt_tree.selectedItems()[0].expt

    db.expt_tree.itemSelectionChanged.connect(sel_change)
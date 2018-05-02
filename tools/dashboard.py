import os, sys, datetime, re, glob, traceback
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
        
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.expt_tree = pg.TreeWidget()
        self.expt_tree.setSortingEnabled(True)
        self.fields = ['timestamp', 'path', 'rig', 'description', 'primary', 'archive', 'backup', 'NAS', 'pipettes.yml', 'site.mosaic', 'DB', 'LIMS', '20x', 'cell map', '63x', 'morphology']
        self.expt_tree.setColumnCount(len(self.fields))
        self.expt_tree.setHeaderLabels(self.fields)
        self.layout.addWidget(self.expt_tree, 0, 0)
        
        self.resize(1000, 900)
        
        self.records = {}
        
        self.poll_thread = PollThread(limit=limit)
        self.poll_thread.update.connect(self.poller_update)
        if no_thread:
            self.poll_thread.poll()  # for local debugging
        else:
            self.poll_thread.start()
        
    def poller_update(self, rec):
        ts = rec['timestamp']
        if ts in self.records:
            item = self.records[ts]['item']
            self.records[ts].update(rec)
        else:
            item = pg.TreeWidgetItem()
            self.expt_tree.addTopLevelItem(item)
            rec['item'] = item
            self.records[ts] = rec
            item.expt = rec['expt']
            
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
                elif val in (False, 'ERROR'):
                    color = fail_color

            item.setText(i, str(val))
            if color is not None:
                item.setBackgroundColor(i, pg.mkColor(color))


class PollThread(QtCore.QThread):
    """Used to check in the background for changes to experiment status.
    """
    update = QtCore.Signal(object)
    
    def __init__(self, limit=0):
        self.limit = limit
        QtCore.QThread.__init__(self)
        
    def run(self):
        while True:
            try:
                self.poll()
            except Exception:
                sys.excepthook(*sys.exc_info())
            break
                
    def poll(self):

        self.session = database.Session()
        expts = {}

        print("loading site paths..")
        #site_paths = glob.glob(os.path.join(config.synphys_data, '*', 'slice_*', 'site_*'))
        
        search_paths = [config.synphys_data]
        for rig_name, rig_path_sets in config.rig_data_paths.items():
            for path_set in rig_path_sets:
                search_paths.extend(list(path_set.values()))

        for path in search_paths:
            print("==============================  Searching %s" % path)
            if not os.path.exists(path):
                print("         skipping; path does not exist.")
                continue
            root_dh = getDirHandle(path)
            for site_dh in self.list_expts(root_dh):
                expt = ExperimentMetadata(path=site_dh.name())
                if expt.timestamp in expts:
                    continue
                expts[expt.timestamp] = expt
                try:
                    rec = self.check(expt)
                    self.update.emit(rec)
                except Exception:
                    rec = {
                        'expt': expt,
                        'path': site_dh.name(),
                        'timestamp': site_dh.info().get('__timestamp__', False),
                        'rig': 'ERROR',
                    }   
                    self.update.emit(rec)
                    traceback.print_exc()
                    print(expt)


                if self.limit > 0 and len(expts) > self.limit:
                    break
            if self.limit > 0 and len(expts) > self.limit:
                break


    def list_expts(self, root_dh):
        for expt in sorted(root_dh.ls(), reverse=True):
            if '@Recycle' in expt:
                continue
            expt_dh = root_dh[expt]
            print(expt_dh.name())
            if not expt_dh.isDir():
                continue
            for slice_name in expt_dh.ls()[::-1]:
                slice_dh = expt_dh[slice_name]
                if not slice_dh.isDir():
                    continue
                print(slice_dh.name())
                for site_name in slice_dh.ls()[::-1]:
                    site_dh = slice_dh[site_name]
                    if not site_dh.isDir():
                        continue
                    yield site_dh

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
                except:
                    description = (org + ' (? genotype)', fail_color)

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
            'backup': False if expt.backup_path is None else (True if os.path.exists(expt.backups_path) else "MISSING"),
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
        self.primary_path = None
        self.archive_path = None
        self.backup_path = None

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
                    # set self.primary_path, self.archive_path, etc.
                    for k,v in path_set.items():
                        setattr(self, k+'_path', os.path.join(v, expt_subpath))
                    break
            if found_paths:
                break

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
        return self.site_info['__timestamp__']

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
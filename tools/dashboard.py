import os, sys, datetime, re, glob
from multipatch_analysis import config, lims
from multipatch_analysis.database import database
from multipatch_analysis.genotypes import Genotype
from acq4.util.DataManager import getDirHandle
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore


fail_color = (255, 200, 200)
pass_color = (200, 255, 200)


class Dashboard(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.expt_tree = pg.TreeWidget()
        self.expt_tree.setSortingEnabled(True)
        self.fields = ['timestamp', 'path', 'rig', 'description', 'primary', 'archive', 'backup', 'NAS', 'pipettes.yml', 'site.mosaic', 'DB', 'LIMS', '20x', 'cell specs', '63x', 'morphology']
        self.expt_tree.setColumnCount(len(self.fields))
        self.expt_tree.setHeaderLabels(self.fields)
        self.layout.addWidget(self.expt_tree, 0, 0)
        
        self.resize(1000, 900)
        
        self.records = {}
        
        self.poll_thread = PollThread()
        self.poll_thread.update.connect(self.poller_update)
        # self.poll_thread.start()
        self.poll_thread.poll()  # for local debugging
        
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
    
    def __init__(self):
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
        
        root_dh = getDirHandle(config.synphys_data)
        
        print(root_dh.name())
        for site_dh in self.list_expts(root_dh):
            expt = ExperimentMetadata(nas_path=site_dh.name())
            if expt.timestamp in expts:
                continue
            expts[expt.timestamp] = expt
            self.check(expt)
            if len(expts) > 50:
                break

    def list_expts(self, root_dh):
        for expt in root_dh.ls()[::-1]:
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
        lims = len(subs) == 1

        rec = {
            'path': expt.site_dh.name(), 
            'timestamp': expt.timestamp, 
            'rig': expt.rig_name, 
            'primary': os.path.exists(expt.rig_path),
            'archive': os.path.exists(expt.archive_path),
            'backup': os.path.exists(expt.backup_path),
            'description': description,
            'pipettes.yml': expt.pipette_file is not None,
            'site.mosaic': expt.mosaic_file is not None,
            'DB': expt.in_database,
            'NAS': expt.nas_path is not None,
            'LIMS': lims,
        }
        
        self.update.emit(rec)


class ExperimentMetadata(object):
    """Handles reading experiment metadata from several possible locations.
    """
    def __init__(self, nas_path=None, rig_path=None, archive_path=None, backup_path=None):
        self._nas_path = None
        self._rig_path = None
        self._archive_path = None
        self._backup_path = None

        if nas_path is not None:
            self._nas_path = nas_path
            self.site_dh = getDirHandle(nas_path)
        elif rig_path is not None:
            self._rig_path = rig_path
            self.site_dh = getDirHandle(rig_path)
        elif archive_path is not None:
            self._archive_path = archive_path
            self.site_dh = getDirHandle(archive_path)
        elif backup_path is not None:
            self._backup_path = backup_path
            self.site_dh = getDirHandle(backup_path)

        self.site_info = self.site_dh.info()
        self._slice_info = None
        self._expt_info = None
        self._specimen_info = None

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
        if self._nas_path is None:
            pass
        return self._nas_path

    @property
    def rig_path(self):
        if self._rig_path is None:
            if self.nas_path is not None:
                self._rig_path = open(os.path.join(self.nas_path, 'sync_source')).read()
        return self._rig_path

    @property
    def archive_path(self):
        if self._archive_path is None:
            if self.nas_path is not None:
                apath = re.sub('mp(\d)', 'mp\1_e', self.rig_path)
                self._archive_path = apath
        return self._archive_path

    @property
    def backup_path(self):
        if self._backup_path is None:
            if self.nas_path is not None:
                self._backup_path = os.path.join('/home/luke/mnt/backup_server', os.path.relpath(self.rig_path, '/home/luke/mnt'))
        return self._backup_path

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
        name = self.expt_info.get('rig_name')
        if name is not None:
            return name

        # infer rig name from paths
        path = None
        if self.archive_path is not None:
            path = self.archive_path
        elif self.rig_path is not None:
            path = self.rig_path
        elif self.backup_path is not None:
            path = self.backup_path
        if path is None:
            return None
        ind = path.lower().index(os.path.sep + 'mp')
        return path[ind+1:ind:4]

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
        return lims.expt_submissions(self.specimen_info['specimen_id'], self.timestamp)


if __name__ == '__main__':
    app = pg.mkQApp()
    pg.dbg()
    db = Dashboard()
    db.show()
    
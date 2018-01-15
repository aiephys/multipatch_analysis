import os, sys, datetime, re, glob
import config
import database
from acq4.util.DataManager import getDirHandle
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore


class Dashboard(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        
        self.layout = QtGui.QGridLayout()
        self.setLayout(self.layout)
        
        self.expt_tree = pg.TreeWidget()
        self.expt_tree.setColumnCount(5)
        self.layout.addWidget(self.expt_tree, 0, 0)
        
        self.resize(1000, 900)
        
        self.records = {}
        
        self.poll_thread = PollThread()
        self.poll_thread.update.connect(self.poller_update)
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
            
        for i, field in enumerate(['timestamp', 'rig', 'has_metadata_qc', 'has_site_mosaic']):
            item.setText(i, str(rec[field]))


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
        print("loading site paths..")
        #site_paths = glob.glob(os.path.join(config.synphys_data, '*', 'slice_*', 'site_*'))
        
        root_dh = getDirHandle(config.synphys_data)
        
        print(root_dh.name())
        for expt in root_dh.ls():
            expt_dh = root_dh[expt]
            print(expt_dh.name())
            if not expt_dh.isDir():
                continue
            for slice_name in expt_dh.ls():
                slice_dh = expt_dh[slice_name]
                if not slice_dh.isDir():
                    continue
                print(slice_dh.name())
                for site_name in slice_dh.ls():
                    site_dh = slice_dh[site_name]
                    if not site_dh.isDir():
                        continue
                    self.check(site_dh)
            
    def check(self, site_dh):
        print("   check %s" % site_dh.name())
        ts = site_dh.info()['__timestamp__']
        date = datetime.datetime.fromtimestamp(ts)
        site_id = date
        
        rig = ''
        #rig = re.search('(MP\d)_', site_dh.name()).groups()[0]
        
        has_meta_qc = 'file_manifest.yml' in site_dh.ls()
        has_site_mosaic = 'site.mosaic' in site_dh.ls()
        
        
        
        #try:
            #expt_entry = database.experiment_from_timestamp(date)
            #expt_steps = expt_entry.submission_data
            #if expt_steps is None:
                #expt_steps = {}
        #except KeyError:
            #expt_entry = None
            #expt_steps = {}
        
        #print("{rig} {date} {uid} {in_db} {in_server}".format(
            #rig=rig,
            #date=date.strftime('%Y-%m-%d'), 
            #uid='%0.2f'%ts,
            #in_db=expt_entry is not None,
            #in_server='raw_data_location' in expt_steps,
        #))
        
        self.update.emit({
            'dh': site_dh.name(), 
            'timestamp': ts, 
            'rig': rig, 
            'has_metadata_qc': has_meta_qc,
            'has_site_mosaic': has_site_mosaic,
        })


if __name__ == '__main__':
    app = pg.mkQApp()
    db = Dashboard()
    db.show()
    
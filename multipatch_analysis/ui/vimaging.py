from acq4.Manager import getManager
from acq4.pyqtgraph.Qt import QtGui, QtCore
import acq4.pyqtgraph as pg
import numpy as np
import scipy.ndimage as ndimage


class VImagingAnalyzer(QtGui.QWidget):
    def __init__(self):
        QtGui.QWidget.__init__(self)
        self.resize(800, 1000)
        
        self.layout = QtGui.QGridLayout()
        self.layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(self.layout)

        self.gw = pg.GraphicsLayoutWidget()
        self.layout.addWidget(self.gw)
        
        self.vb1 = self.gw.addViewBox()
        self.img1 = pg.ImageItem()
        self.vb1.addItem(self.img1)
        
        self.vb2 = self.gw.addViewBox(row=1, col=0)
        self.img2 = pg.ImageItem()
        self.vb2.addItem(self.img2)
        
        for vb in (self.vb1, self.vb2):
            vb.invertY()
            vb.setAspectLocked(True)

        self.plt1 = self.gw.addPlot(row=2, col=0)
        self.plt1.setLabels(bottom=('time', 's'))
        self.plt1_items = []
        
        self.rois = [pg.RectROI([0, 0], [10, 10], pen=(i,4)) for i in range(2)]
        for roi in self.rois:
            self.vb1.addItem(roi)
            roi.sigRegionChanged.connect(self.roi_changed)
            roi.sigRegionChangeFinished.connect(self.update_sequence_analysis)
            
        self.plt2 = self.gw.addPlot(row=3, col=0)
        self.plt2.setXLink(self.plt1)
        self.plt2.setLabels(bottom=('time', 's'))
            
        self.base_time_rgn = pg.LinearRegionItem([0.0, 0.1])
        self.test_time_rgn = pg.LinearRegionItem([0.1, 0.11])
        for time_rgn in (self.base_time_rgn, self.test_time_rgn):
            self.plt1.addItem(time_rgn)
            time_rgn.sigRegionChanged.connect(self.time_rgn_changed)
            time_rgn.sigRegionChangeFinished.connect(self.update_sequence_analysis)
        self.clamp_plots = []
        
        self.plt3 = self.gw.addPlot(row=4, col=0)
        self.plt3.setLabels(left="dF / F", bottom="trial")
        
        self.show()

    def load_data(self, seqDir):
        with pg.BusyCursor():
            self._load_data(seqDir)

    def _load_data(self, seqDir):
        man = getManager()
        model = man.dataModel

        # read all image data
        self.img_data = model.buildSequenceArray(
            seqDir, 
            lambda dh: dh['Camera']['frames.ma'].read(),
            join=False).asarray()
        seqParams = list(model.listSequenceParams(seqDir).items())
        if self.img_data.ndim == 1:
            self.img_data = self.img_data[np.newaxis, :]
            seqParams.insert(0, (None, [0]))

        transpose = seqParams[0][0] == ('protocol', 'repetitions')
        if transpose:
            self.img_data = np.swapaxes(self.img_data, 0, 1)
            seqParams = seqParams[::-1]
        self.seqParams = seqParams

        if seqParams[0][0] is None:
            self.seqColors = [pg.mkColor('w')]
        else:
            nSeq = len(seqParams[0][1])
            self.seqColors = [pg.intColor(i, nSeq*1.6) for i in range(nSeq)]

        # cull out truncated recordings :(
        self.img_data = [[d for d in row if d.xvals('Time')[-1] > 150e-3] for row in self.img_data]

        # crop / concatenate
        img_len = min([min([d.shape[0] for d in row]) for row in self.img_data])
        self.img_data = [[d[:img_len] for d in row] for row in self.img_data]
        self.img_arrays = [np.concatenate([d.asarray()[np.newaxis, ...] for d in row], axis=0) for row in self.img_data]

        # average
        self.img_mean = [img_arr.mean(axis=0) for img_arr in self.img_arrays]

        for p in self.clamp_plots:
            self.plt2.removeItem(p)
            
        # read all clamp data
        first_subdir = seqDir[seqDir.ls()[0]]
        clamp_file = first_subdir['Clamp2.ma']
        if clamp_file.exists():
            self.clamp_mode = model.getClampMode(clamp_file)
            chan = 'command' if self.clamp_mode == 'VC' else 'primary'
            self.clamp_data = model.buildSequenceArray(
                seqDir, 
                lambda dh: dh['Clamp2.ma'].read()['Channel': chan],
                join=False).asarray()
            if self.clamp_data.ndim == 1:
                self.clamp_data = self.clamp_data[np.newaxis, :]
            if transpose:
                self.clamp_data = np.swapaxes(self.clamp_data, 0, 1)

            self.plt2.setLabels(left=('Vm', 'V'))
            
            for i in range(self.clamp_data.shape[0]):
                for j in range(self.clamp_data.shape[1]):
                    trace = self.clamp_data[i,j]
                    pen = self.seqColors[i]
                    p = self.plt2.plot(trace.xvals('Time'), trace.asarray(), antialias=True, pen=pen)
                    self.clamp_plots.append(p)
            self.plt2.show()
        else:
            self.clamp_mode = None
            self.plt2.hide()

        self.img_t = self.img_data[0][0].xvals('Time')

        self.img1.setImage(self.img_mean[-1].mean(axis=0))

        self.roi_changed(None)
        self.time_rgn_changed(None)
        self.update_sequence_analysis()

    def roi_changed(self, roi):
        for item in self.plt1_items:
            self.plt1.removeItem(item)

        roi1, roi2 = self.rois
        
        for i, img_data in enumerate(self.img_arrays):
            color = self.seqColors[i]
            color2 = pg.mkColor(color)
            color2.setAlpha(40)

            rgn1 = roi1.getArrayRegion(img_data, self.img1, axes=(2, 3)).mean(axis=2).mean(axis=2)
            rgn2 = roi2.getArrayRegion(img_data, self.img1, axes=(2, 3)).mean(axis=2).mean(axis=2)
            dif = rgn1 - rgn2
            
            difmean = dif.mean(axis=0)
            baseline = np.median(difmean[:10])
            # plot individual examples only for last parameter
            if i == len(self.img_arrays)-1:
                for j in range(dif.shape[0]):
                    offset = baseline - np.median(dif[j, :10])
                    self.plt1_items.append(self.plt1.plot(self.img_t, dif[j] + offset, pen=color2, antialias=True))
            # plot average
            self.plt1_items.append(self.plt1.plot(self.img_t, difmean, pen=color, antialias=True))
            self.plt1_items[-1].setZValue(10)
        
    def time_rgn_changed(self, rgn):
        img_data = self.img_arrays[-1]
        img_t = self.img_data[-1][0].xvals('Time')
        
        base_starti, base_stopi, test_starti, test_stopi = self.time_indices(img_t)

        base = img_data[:, base_starti:base_stopi].mean(axis=0).mean(axis=0)
        test = img_data[:, test_starti:test_stopi].mean(axis=0).mean(axis=0)

        dff = ndimage.median_filter(test - base, 10)
        self.img2.setImage(dff)

    def time_indices(self, time_vals):
        base_start, base_stop = self.base_time_rgn.getRegion()
        base_starti = np.argwhere(time_vals >= base_start)[0,0]
        base_stopi = np.argwhere(time_vals >= base_stop)[0,0]
        
        test_start, test_stop = self.test_time_rgn.getRegion()
        test_starti = np.argwhere(time_vals >= test_start)[0,0]
        test_stopi = np.argwhere(time_vals >= test_stop)[0,0]
        
        return base_starti, base_stopi, test_starti, test_stopi

    def update_sequence_analysis(self):
        xvals = []
        yvals = []
        brushes = []
        avg_x = []
        avg_y = []

        for i in range(len(self.img_data)):

            img_data = self.img_arrays[i]
            img_t = self.img_data[i][0].xvals('Time')
            base_starti, base_stopi, test_starti, test_stopi = self.time_indices(img_t)
            
            base = img_data[:, base_starti:base_stopi]
            base_mean = base.mean(axis=1)
            test = img_data[:, test_starti:test_stopi]
            test_mean = test.mean(axis=1)
            
            roi1, roi2 = self.rois
            base_rgn1 = roi1.getArrayRegion(base_mean, self.img1, axes=(1, 2))
            base_rgn2 = roi2.getArrayRegion(base_mean, self.img1, axes=(1, 2))
            test_rgn1 = roi1.getArrayRegion(test_mean, self.img1, axes=(1, 2))

            background = base_rgn2.mean(axis=1).mean(axis=1)
            baseline = base_rgn1.mean(axis=1).mean(axis=1)
            signal = test_rgn1.mean(axis=1).mean(axis=1)
            dff = (signal - baseline) / (baseline - background)

            x = self.seqParams[0][1][i]
            xvals.extend([x] * len(dff))
            yvals.extend(list(dff))
            color = pg.mkColor(self.seqColors[i])
            color.setAlpha(50)
            brushes.extend([pg.mkBrush(color)] * len(dff))
            avg_x.append(x)
            avg_y.append(dff.mean())

        self.plt3.clear()
        if self.seqParams[0][0] is None:
            self.plt3.plot(range(len(xvals)), yvals, pen='w')
        else:
            self.plt3.plot(xvals, yvals, symbol='o', pen=None, symbolBrush=brushes, symbolPen=None)
            self.plt3.plot(avg_x, avg_y, pen='w', antialias=True)

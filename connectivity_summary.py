# *-* coding: utf-8 *-*
"""
I have a summary of ALM multipatch experiments stored on workflowy and exported to expt_summary.txt.

This file describes, for each experiment:
   - which cre labels were used
   - which cells are +/- for which labels
   - which cells had usable recordings
   - which cells were connected

The purpose of this script is to parse the exported file and return a summary
of total connections detected / probed for each cell type pair.

"""
from __future__ import print_function, division
import os, re, sys, traceback, glob, json, warnings, datetime, argparse, pickle
import numpy as np
import scipy.optimize, scipy.stats
import pyqtgraph as pg
import pyqtgraph.configfile
import allensdk_internal.core.lims_utilities as lims


ALL_CRE_TYPES = ['sst', 'pvalb', 'tlx3', 'sim1', 'rorb', 'vip', 'ntsr1', 'chrna2', 'rbp4', 'chat', 'ctgf', 'ndnf', 'glt25d2', 'htr3a', 'nos1', ]
ALL_LABELS = ['biocytin', 'af488', 'cascade_blue']


def arg_to_date(arg):
    if arg is None:
        return None
    parts = re.split('\D+', arg)
    return datetime.date(*map(int, parts))




# read exported file

#steph = open('/home/luke/mnt/mp1/Steph/DataSummary', 'r').readlines()
#pasha = open('/home/luke/mnt/mp2/data/Pasha/connection_analysis', 'r').readlines()
#alex = open('/home/luke/mnt/mp3/data/Alex/connection analysis', 'r').readlines()
#lines = steph + pasha + alex



# first parse indentation to generate a hierarchy of strings
# Note: summary was exported in plain text mode, which can be difficult to parse.
# Might want to switch to XML if this becomes an issue.

class Entry(object):
    def __init__(self, line, parent, file, lineno):
        self.indentation = indentation(line)
        self.lines = []
        self.add_line(line)
        self.parent = parent
        self.file = file
        self.lineno = lineno
        self.children = []
        if parent is not None:
            parent.add_child(self)
        
    def add_line(self, line):
        if indentation(line) != self.indentation:
            raise IndentationError(line)
        self.lines.append(line[self.indentation:].rstrip())
        
    def add_child(self, child):
        self.children.append(child)

    def print_tree(self):
        print("\n".join([('    '*self.indentation)+l for l in self.lines]))
        for ch in self.children:
            ch.print_tree()

        
def indentation(line):
    return len(line) - len(line.lstrip('- '))



# Now go through each experiment and read cell type / connectivity data

class Experiment(object):
    def __init__(self, entry):
        self.entry = entry
        self.connections = []
        self._region = None
        self._summary = None
        self._view = None
        self._slice_info = None
        self._lims_record = None
        self._site_path = None
        self._probed = None

        try:
            self.expt_id = entry.lines[0]
            self.cells = {x:Cell(self,x) for x in range(1,9)}
        except Exception as exc:
            Exception("Error parsing experiment: %s\n%s" % (self, exc.args[0]))
        
        for ch in entry.children:
            try:
                if ch.lines[0] == 'Labeling':
                    self.parse_labeling(ch)
                elif ch.lines[0] == 'Cell QC':
                    self.parse_qc(ch)
                elif ch.lines[0] == 'Connections':
                    self.parse_connections(ch)
                elif ch.lines[0] == 'Conditions':
                    continue
                elif ch.lines[0].startswith('Region '):
                    self._region = ch.lines[0][7:]
                elif ch.lines[0].startswith('Site path '):
                    p = os.path.abspath(os.path.join(os.path.dirname(self.entry.file), ch.lines[0][10:]))
                    if not os.path.isdir(p):
                        raise Exception("Invalid site path: %s" % p)
                    self._site_path = p
                else:
                    raise Exception('Invalid experiment entry "%s"' % ch.lines[0])
                
            except Exception as exc:
                raise Exception("Error parsing %s for experiment: %s\n%s" % (ch.lines[0], self, exc.args[0]))
            
        # gather lists of all labels and cre types
        cre_types = set()
        labels = set()
        for cell in self.cells.values():
            cre_types |= set(cell.labels.keys()) & set(ALL_CRE_TYPES)
            labels |= set(cell.labels.keys()) & set(ALL_LABELS)
        self.cre_types = sorted(list(cre_types), key=lambda x: ALL_CRE_TYPES.index(x))
        self.labels = sorted(list(labels), key=lambda x: ALL_LABELS.index(x))
        
        # make sure all cells have information for all labels
        for cell in self.cells.values():
            for label in labels:
                assert label in cell.labels
            for cre in cre_types:
                assert cre in cell.labels
                
        # read cell positions from mosaic files
        try:
            self.load_cell_positions()
        except Exception as exc:
            print("Warning: Could not load cell positions for %s:\n    %s" % (self, exc.args[0]))
                
        # pull donor/specimen info from LIMS
        self.age
        # check for a single NWB file
        self.nwb_file

    def parse_labeling(self, entry):
        """
        "Labeling" section should look like:
        
            Labeling:
                sim1: 1+ 2- 3x 4x+ 5+? 6?
                biocytin: ...
                af488: 1+ 2+ 3x 4- 5? 6+
                cascade_blue: ...
                
        This example has the following interpretation:
        
            1+   Cell 1 is reporter-positive and dye filled
            2-   Cell 2 is reporter-negative and dye filled   
            3x   Cell 3 type cannot be determined (no pipette tip found)
            4x+  Cell 4 was not dye filled, but pipette tip appears to be touching cre-positive cell
            5+?  Cell 5 looks like cre-positive, but image is somewhat unclear
            6?   Cell 6 is dye filled, but cre type is ambiguous
        """
        for ch in entry.children:
            line = ch.lines[0]
            # line looks like "sim1: 1+ 2-'
            
            parts = re.split('\s+', line)
            
            # first part is label / cre type and a colon
            assert parts[0].endswith(':')
            cre = parts[0][:-1]
            if not (cre in ALL_LABELS or cre in ALL_CRE_TYPES):
                raise Exception("Invalid label or cre type: %s" % cre)
            
            # parse the remainder of the line
            if len(parts[1:]) == 1 and parts[1].strip() == '?':
                # no data
                continue

            for part in parts[1:]:
                m = re.match('(\d+)(x)?(\+|\-)?(\?)?', part)
                if m is None:
                    raise Exception('invalid label record: %s' % part)
                grps = m.groups()
                cell_id = int(grps[0])
                cell = self.cells[cell_id]
                absent = grps[1] == 'x'
                positive = grps[2]
                uncertain = grps[3] == '?'
                
                assert cre not in cell.labels
                cell.labels[cre] = positive
        
    def parse_qc(self, entry):
        """Parse cell quality control. Looks like:
        
            Holding: 1- 2+ 3- 4- 6/ 7+
            Access: 1- 2/ 3+ 4- 6/ 7/
            Spiking: 1- 2+ 3+ 4- 6+ 7+
        
        Where + means pass, / means borderline pass, - means fail, ? means unknown
        """
        for ch in entry.children:
            parts = re.split('\s+', ch.lines[0].strip())
            for part in parts[1:]:
                m = re.match(r'(\d+)(\+|\/|\-|\?)', part)
                if m is None:
                    raise Exception('Invalid cell QC string "%s"' % part)
                cell_id = int(m.groups()[0])
                val = m.groups()[1]
                cell = self.cells[cell_id]
                
                if parts[0] == 'Holding:':
                    cell.holding_qc = val in '+/'
                elif parts[0] == 'Access:':
                    cell.access_qc = val in '+/'                    
                elif parts[0] == 'Spiking:':
                    cell.spiking_qc = val in '+/'                    
                else:
                    raise Exception("Invalid Cell QC line: %s" % ch.lines[0])
    
    def parse_connections(self, entry):
        if len(entry.children) == 0 or entry.children[0].lines[0] == 'None':
            return
        for con in entry.children:
            m = re.match(r'(\d+)\s*->\s*(\d+)\s*(\??)\s*(.*)', con.lines[0].strip())
            if m is None:
                raise Exception("Invalid connection: %s" % con.lines[0])
            
            if m.groups()[2] == '?':
                # ignore questionable connections
                continue
            self.connections.append((int(m.groups()[0]), int(m.groups()[1])))

    def summary(self):
        """Return a structure summarizing (non)connectivity in the experiment.
        
        Looks like:
        
            {(pre_type, post_type): [connected, unconnected, cdist, udist], ...}
        """
        if self._summary is None:
            csum = {}
            for i, j in self.connections_probed:
                ci, cj = self.cells[i], self.cells[j]
                typ = (ci.cre_type, cj.cre_type)
                if typ not in csum:
                    csum[typ] = [0, 0, [], []]
                if (i, j) in self.connections:
                    csum[typ][0] += 1
                    csum[typ][2].append(ci.distance(cj))
                else:
                    csum[typ][1] += 1
                    csum[typ][3].append(ci.distance(cj))
            self._summary = csum
        return self._summary
    
    @property
    def region(self):
        return 'V1' if (not hasattr(self, '_region') or self._region is None) else self._region
    
    @property
    def connections_probed(self):
        """A list of probed connections (pre_cell, post_cell) that passed QC.
        """
        if self._probed is None:
            probed = []
            for i,ci in self.cells.items():
                for j,cj in self.cells.items():
                    if i == j:
                        continue
                    if ci.spiking_qc is not True:
                        # presynaptic cell failed spike QC; ignore
                        continue
                    if cj.pass_qc is not True:
                        # postsynaptic cell failed ephys qc; ignore
                        continue
                    if ci.cre_type is None or cj.cre_type is None:
                        # indeterminate cell types; ignore
                        #print("Ignore probe (ind. cell type) %s:%d-%d" % (self.expt_id, i, j))
                        #if (i, j) in self.connections:
                            #print("    --> connected!")
                        continue
                    probed.append((i, j))
            self._probed = probed
        return self._probed
    
    @property
    def n_connections_probed(self):
        return len(self.connections_probed)

    @property
    def n_connections(self):
        return sum([x[0] for x in self.summary().values()])
        
    def load_cell_positions(self):
        """Load cell positions from external file.
        """
        sitefile = self.mosaic_file
        mosaic = json.load(open(sitefile))
        marker_items = [i for i in mosaic['items'] if i['type'] == 'MarkersCanvasItem']
        if len(marker_items) == 0:
            raise TypeError("No cell markers found in site mosaic file %s" % sitefile)
        elif len(marker_items) > 1:
            raise TypeError("Multiple marker items found in site mosaic file %s" % sitefile)
        cells = marker_items[0]['markers']
        for name, pos in cells:
            m = re.match("\D+(\d+)", name)
            cid = int(m.group(1))
            self.cells[cid].position = pos

    @property
    def mosaic_file(self):
        """Path to site mosaic file
        """
        sitefile = os.path.join(self.path, "site.mosaic")
        if not os.path.isfile(sitefile):
            sitefile = os.path.join(os.path.split(self.path)[0], "site.mosaic")
        if not os.path.isfile(sitefile):
            mosaicfiles = [f for f in os.listdir(self.path) if f.endswith('.mosaic')]
            if len(mosaicfiles) == 1:
                sitefile = os.path.join(self.path, mosaicfiles[0])
        if not os.path.isfile(sitefile):
            # print(os.listdir(self.path))
            # print(os.listdir(os.path.split(self.path)[0]))
            raise Exception("No site mosaic found for %s" % self)
        return sitefile
        
    @property
    def path(self):
        """Filesystem path to the root of this experiment.
        """
        if self._site_path is None:
            date, slice, site = self.expt_id.split('-')
            root = os.path.dirname(self.entry.file)
            paths = [
                os.path.join(root, date+"_000", "slice_%03d"%int(slice), "site_%03d"%int(site)),
                os.path.join(root, 'V1', date+"_000", "slice_%03d"%int(slice), "site_%03d"%int(site)),
                os.path.join(root, 'ALM', date+"_000", "slice_%03d"%int(slice), "site_%03d"%int(site)),
            ]
            for path in paths:
                if os.path.isdir(path):
                    self._site_path = path
                    break
            if self._site_path is None:
                raise Exception("Cannot find filesystem path for experiment %s. Attempted paths:\n%s" % (self, "\n".join(paths)))
        return self._site_path
    
    def __repr__(self):
        return "<Experiment %s (%s:%d)>" % (self.expt_id, self.entry.file, self.entry.lineno)

    @property
    def slice_info(self):
        if self._slice_info is None:
            index = os.path.join(os.path.split(self.path)[0], '.index')
            if not os.path.isfile(index):
                raise TypeError("Cannot find index file (%s) for experiment %s" % (index, self))
            self._slice_info = pg.configfile.readConfigFile(index)['.']
        return self._slice_info

    @property
    def nwb_file(self):
        p = self.path
        files = glob.glob(os.path.join(p, '*.nwb'))
        if len(files) == 0:
            files = glob.glob(os.path.join(p, '*.NWB'))
        if len(files) == 0:
            raise Exception("No NWB file found for %s" % self)
        elif len(files) > 1:
            raise Exception("Multiple NWB files found for %s" % self)
        return files[0]

    @property
    def specimen_id(self):
        return self.slice_info['specimen_ID'].strip()

    @property
    def age(self):
        age = self.lims_record.get('days', 0)
        if age == 0:
            raise Exception("Donor age not set in LIMS for specimen %s" % self.specimen_id)
            # data not entered in to lims
            age = (self.date - self.birth_date).days
        return age

    @property
    def birth_date(self):
        bd = self.lims_record['date_of_birth']
        return datetime.date(bd.year, bd.month, bd.day)

    @property
    def lims_record(self):
        if self._lims_record is None:
            sid = self.specimen_id
            q = """
                select d.date_of_birth, ages.days from donors d
                join specimens sp on sp.donor_id = d.id
                join ages on ages.id = d.age_id
                where sp.name  = '%s'
                limit 2""" % sid
            r = lims.query(q)
            if len(r) != 1:
                raise Exception("LIMS lookup for specimen %s returned %d results" % (sid, len(r)))
            self._lims_record = r[0]
        return self._lims_record

    @property
    def date(self):
        y,m,d = self.expt_id.split('-')[0].split('.')
        return datetime.date(int(y), int(m), int(d))

    def show(self):
        if self._view is None:
            pg.mkQApp()
            self._view_widget = pg.GraphicsLayoutWidget()
            self._view = self._view_widget.addViewBox(0, 0)
            v = self._view
            cell_ids = sorted(self.cells.keys())
            pos = np.array([self.cells[i].position[:2] for i in cell_ids])
            if len(self.connections) == 0:
                adj = np.empty((0,2), dtype='int')
            else:
                adj = np.array(self.connections) - 1
            colors = []
            for cid in cell_ids:
                cell = self.cells[cid]
                color = [0, 0, 0]
                for i,cre in enumerate(self.cre_types):
                    if cell.labels[cre] == '+':
                        color[i] = 255
                colors.append(color)
            brushes = [pg.mkBrush(c) for c in colors]
            print(pos)
            print(adj)
            print(colors)
            self._graph = pg.GraphItem(pos=pos, adj=adj, size=30, symbolBrush=brushes)
            v.addItem(self._graph)
        self._view_widget.show()
        
        
class Cell(object):
    def __init__(self, expt, cell_id):
        self.expt = expt
        self.cell_id = cell_id
        self.access_qc = None
        self.holding_qc = None
        self.spiking_qc = None
        self.labels = {}
        self.position = None
    
    @property
    def pass_qc(self):
        """True if cell passes QC.
        """
        if self.access_qc is True and self.holding_qc is True:
            return True
        elif self.access_qc is False or self.holding_qc is False:
            return False
        
        # None means cell is not present in ephys data
        return None

    @property
    def cre_type(self):
        """Cre type string for this cell.
        
        If the cell is reporter-negative then cre_type is 'unk'.
        If the cell has ambiguous or missing data then cre_type is None.
        """
        default = 'unknown'
        ct = None
        for label,pos in self.labels.items():
            if label in ALL_LABELS:
                continue
            if pos == '+':
                if ct not in (None, default):
                    raise Exception("%s has multiple labels!" % self)
                ct = label
            elif pos == '-':
                if ct is not None:
                    continue
                ct = default
        return ct

    @property
    def label_type(self):
        """fluorescent type string for this cell.
        
        If the cell is reporter-negative then cre_type is 'unk'.
        If the cell has ambiguous or missing data then cre_type is None.
        """
        default = 'unknown'
        ct = None
        for label,pos in self.labels.items():
            if label in ALL_CRE_TYPES or label == 'biocytin':
                continue
            if pos == '+':
                if ct not in (None, default):
                    raise Exception("%s has multiple labels!" % self)
                ct = label
            elif pos == '-':
                if ct is not None:
                    continue
                ct = default
        return ct

    def distance(self, cell):
        """Return distance between cells, or nan if positions are not defined.
        """
        p1 = self.position
        p2 = cell.position
        if p1 is None or p2 is None:
            return np.nan
        return ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)**0.5

    def __repr__(self):
        return "<Cell %s:%d>" % (self.expt.expt_id, self.cell_id)


class ExperimentList(object):
    def __init__(self, expts=None):
        if expts is None:
            expts = []
        self._expts = expts
        self.start_skip = []
        self.stop_skip = []

    def load(self, filename):
        if filename.endswith('.pkl'):
            self._load_pickle(filename)
        else:
            self._load_text(filename)

    def _load_pickle(self, filename):
        el = pickle.load(open(filename))
        self._expts.extend(el._expts)
        self.sort()
    
    def _load_text(self, filename):
        root = Entry('', None, None, None)
        root.indentation = -1
        current = root

        lines = open(filename, 'r').readlines()

        for i,line in enumerate(lines):
            if line.lstrip().startswith('#'):
                continue
            
            if line.strip() == '':
                continue
            ind = indentation(line)
            if ind > current.indentation:
                ch = Entry(line, parent=current, file=f, lineno=i)
                current = ch
            else:
                while current.indentation > ind:
                    current = current.parent
                ch = Entry(line, parent=current.parent, file=f, lineno=i)
                current = ch
                continue
            
        # Parse experiment data
        errs = []
        for entry in root.children:
            try:
                expt = Experiment(entry)
            except Exception as exc:
                errs.append((entry, sys.exc_info()))
                continue
            
            self._expts.append(expt)

        if len(errs) > 0:
            print("Errors loading %d experiments:" % len(errs))
            for entry, exc in errs:
                print("=======================")
                print("\n".join(entry.lines))
                traceback.print_exception(*exc)
                print("")

        self.sort()

    def select(self, start=None, stop=None, region=None):
        expts = []
        start_skip = []
        stop_skip = []
        for ex in self._expts:
            # filter experiments by start/stop dates
            if args.start is not None and ex.date < args.start:
                continue
            elif args.stop is not None and ex.date > args.stop:
                continue
            elif args.region is not None and ex.region != args.region:
                continue
            else:
                expts.append(ex)
    
        el = ExperimentList(expts)
        el.start_skip = self.start_skip + start_skip
        el.stop_skip = stop_skip + self.stop_skip
        return el
    
    def __getitem__(self, item):
        return self._expts[item]
    
    def __len__(self):
        return len(self._expts)
    
    def __iter__(self):
        return self._expts.__iter__()
    
    def append(self, expt):
        self._expts.append(expt)

    def sort(self, key=lambda expt: expt.expt_id, **kwds):
        self._expts.sort(key=key, **kwds)

    def check(self):
        # sanity check: all experiments should have cre and fl labels
        for expt in expts:
            # make sure we have at least one non-biocytin label and one cre label
            if len(expt.cre_types) < 1:
                print("Warning: Experiment %s has no cre-type labels" % expt.expt_id)
            if len(expt.labels) < 1 or expt.labels == ['biocytin']:
                print("Warning: Experiment %s has no fluorescent labels" % expt.expt_id)
            if expt.region is None:
                print("Warning: Experiment %s has no region" % expt.expt_id)
    
    def distance_plot(self, pre_type, post_type, plot=None, color=(100, 100, 255)):
        # get all connected and unconnected distances for pre->post
        probed = []
        connected = []
        for expt in self:
            for i,j in expt.connections_probed:
                ci, cj = expt.cells[i], expt.cells[j]
                if ci.cre_type != pre_type or cj.cre_type != post_type:
                    continue
                dist = ci.distance(cj)
                probed.append(dist)
                connected.append((i, j) in expt.connections)
        connected = np.array(connected).astype(float)
        probed = np.array(probed)
        pts = np.vstack([probed, connected]).T

        # scatter points a bit
        conn = pts[:,1] == 1
        unconn = pts[:,1] == 0
        if np.any(conn):
            cscat = pg.pseudoScatter(pts[:,0][conn], spacing=10e-6, bidir=False)
            pts[:,1][conn] -= cscat * 0.2 / cscat.max()
        if np.any(unconn):
            uscat = pg.pseudoScatter(pts[:,0][unconn], spacing=10e-6, bidir=False)
            pts[:,1][unconn] -= uscat * 0.2 / uscat.max()

        # scatter plot connections probed
        if plot is None:
            plot = pg.plot()
            
        plot.setLabels(bottom=('distance', 'm'), left='connection probability')
            
        color2 = color + (100,)
        scatter = plot.plot(pts[:,0], pts[:,1], pen=None, symbol='o', labels={'bottom': ('distance', 'm')}, symbolBrush=color2, symbolPen=None)

        # use a sliding window to plot the proportion of connections found along with a 95% confidence interval
        # for connection probability
        def binomial_ci(p, n, alpha=0.05 ):
            """
            Two-sided confidence interval for a binomial test.

            If after n trials we obtain p successes, find c such that

            P(k/n < p/n; theta = c) = alpha

            where k/N is the proportion of successes in the set of trials,
            and theta is the success probability for each trial. 
            
            Source: http://stackoverflow.com/questions/13059011/is-there-any-python-function-library-for-calculate-binomial-confidence-intervals
            """
            upper_fn = lambda c: scipy.stats.binom.cdf(p, n, c) - alpha
            lower_fn = lambda c: scipy.stats.binom.cdf(p, n, c) - (1.0 - alpha)
            return scipy.optimize.bisect(lower_fn, 0, 1), scipy.optimize.bisect(upper_fn, 0, 1)
        
        window = 40e-6
        spacing = window / 4.0
        xvals = np.arange(window / 2.0, 500e-6, spacing)
        upper = []
        lower = []
        prop = []
        ci_xvals = []
        for x in xvals:
            minx = x - window / 2.0
            maxx = x + window / 2.0
            # select points inside this window
            mask = (probed >= minx) & (probed <= maxx)
            pts_in_window = connected[mask]
            # compute stats for window
            n_probed = pts_in_window.shape[0]
            n_conn = pts_in_window.sum()
            if n_probed == 0:
                prop.append(np.nan)
            else:
                prop.append(n_conn / n_probed)
                ci = binomial_ci(n_conn, n_probed)
                lower.append(ci[0])
                upper.append(ci[1])
                ci_xvals.append(x)
        
        # plot connection probability and confidence intervals
        color2 = [c / 3.0 for c in color]
        mid_curve = plot.plot(xvals, prop, pen=color, antialias=True)
        upper_curve = plot.plot(ci_xvals, upper, pen=color2, antialias=True)
        lower_curve = plot.plot(ci_xvals, lower, pen=color2, antialias=True)
        upper_curve.setVisible(False)
        lower_curve.setVisible(False)
        color2 = color + (50,)
        fill = pg.FillBetweenItem(upper_curve, lower_curve, brush=color2)
        fill.setZValue(-10)
        plot.addItem(fill, ignoreBounds=True)
        
        return scatter, mid_curve, lower_curve, upper_curve, fill 

    def matrix(self, rows, cols, size=50):
        w = pg.GraphicsLayoutWidget()
        v = w.addViewBox()
        v.setAspectLocked()
        v.invertY()
        
        colormap = pg.ColorMap(
            [0, 0.01, 0.03, 0.1, 0.3, 1.0],
            [(0,0,100), (80,0,80), (140,0,0), (255,100,0), (255,255,100), (255,255,255)],
        )
        default = (0, 0, 0)
        
        summary = self.connectivity_summary()
        
        w.matrix_items = []
        
        for i,row in enumerate(rows):
            txt = pg.QtGui.QGraphicsTextItem(row)
            br = txt.boundingRect()
            txt.setPos(-br.width() - 10, i * size + size/2. - br.center().y())
            txt.setDefaultTextColor(pg.mkColor('w'))
            v.addItem(txt)
            w.matrix_items.append(txt)
        
        for i,col in enumerate(cols):
            txt = pg.QtGui.QGraphicsTextItem(col)
            br = txt.boundingRect()
            txt.setPos(i * size + size/2 - br.center().x(), -br.height() - 10)
            txt.setDefaultTextColor(pg.mkColor('w'))
            v.addItem(txt)
            w.matrix_items.append(txt)
        
        for i,row in enumerate(rows):
            for j,col in enumerate(cols):
                if (row, col) in summary:
                    conn, uconn, cdist, udist = summary[(row, col)]
                    color = colormap.map(conn/(conn + uconn))
                else:
                    conn, uconn = 0, 0
                    color = default
                x = j * size
                y = i * size
                rect = pg.QtGui.QGraphicsRectItem(x, y, size, size)
                rect.setBrush(pg.mkBrush(color))
                rect.setZValue(-10)
                w.matrix_items.append(rect)
                v.addItem(rect)
                
                txt = pg.QtGui.QGraphicsTextItem("%d/%d" % (conn, conn+uconn))
                br = txt.boundingRect()
                txt.setPos(x + size/2 - br.center().x(), y + size/2 - br.center().y())
                c = 'w' if sum(color) < 300 else 'k'
                if conn == uconn == 0:
                    c = 0.3
                txt.setDefaultTextColor(pg.mkColor(c))
                v.addItem(txt)
                w.matrix_items.append(txt)
                
        # colormap is logarithmic; remap to linear for legend
        colors = colormap.color
        x = np.linspace(0, 1, len(colors))
        cmap2 = pg.ColorMap(x, colors)
        legend = pg.GradientLegend([25, 300], [-20, -30])
        legend.setGradient(cmap2.getGradient())
        legend.setLabels({'%d'%int(a*100):b for a,b in zip(colormap.pos, x)})
        v.addItem(legend)
        w.show()
        self.matrix_widget = w

    def n_connections_probed(self):
        """Return (total_probed, total_connected) for all experiments in this list.
        """
        tot_probed = 0
        tot_connected = 0
        for expt in self:
            tot_probed += expt.n_connections_probed
            tot_connected += expt.n_connections
        return tot_probed, tot_connected
        
    def print_expt_summary(self, list_stims=False):
        fields = ['# probed', '# connected', 'age', 'cre types']
        if list_stims:
            fields.append('stim sets')
        print("----------------------------------------------------------")
        print("  Experiment Summary  (%s)" % ', '.join(fields))
        print("----------------------------------------------------------")

        if len(self.start_skip) > 0:
            print("[ skipped %d earlier experiments ]" % len(self.start_skip))
        tot_probed = 0
        tot_connected = 0
        ages = []
        for i,expt in enumerate(self):
            n_p = expt.n_connections_probed
            n_c = expt.n_connections
            tot_probed += n_p
            tot_connected += n_c
            ages.append(expt.age)
            
            fmt = "%d: %s:  \t%d\t%d\t%d\t%s"
            fmt_args = [i, expt.expt_id, n_p, n_c, expt.age, ', '.join(expt.cre_types)]
            
            # get list of stimuli
            if list_stims:
                from neuroanalysis.miesnwb import MiesNwb
                nwb = MiesNwb(expt.nwb_file)
                stims = []
                for srec in nwb.contents:
                    stim = srec.recordings[0].meta['stim_name']
                    if stim.startswith('PulseTrain_'):
                        stim = stim[11:]
                    if stim.endswith('_DA_0'):
                        stim = stim[:-5]
                    if stim not in stims:
                        stims.append(stim)
                nwb.close()
                
                # sort by frequency
                def freq(stim):
                    m = re.match('(.*)(\d+)Hz', stim)
                    if m is None:
                        return (0, 0)
                    else:
                        return (len(m.groups()[0]), int(m.groups()[1]))
                stims.sort(key=freq)
                
                fmt += "\t%s"
                fmt_args.append(', '.join(stims))
            
            print(fmt % tuple(fmt_args))

        if len(self.stop_skip) > 0:
            print("[ skipped %d later experiments ]" % len(self.stop_skip))
        print("")

        print("Mean age: %0.1f" % np.mean(ages))
        print("")

    def connectivity_summary(self):
        summary = {}
        for expt in self:
            for k,v in expt.summary().items():
                if k not in summary:
                    summary[k] = [0, 0, [], []]
                summary[k][0] += v[0]
                summary[k][1] += v[1]
                summary[k][2].extend(v[2])
                summary[k][3].extend(v[3])
        return summary
        
    def print_connectivity_summary(self):
        print("-------------------------------------------------------------")
        print("     Connectivity  (# connected/probed, % connectivity, cdist, udist, adist)")
        print("-------------------------------------------------------------")
        
        tot_probed, tot_connected = self.n_connections_probed()
        
        summary = self.connectivity_summary()

        with warnings.catch_warnings():  # we expect warnings when nanmean is called on an empty list
            warnings.simplefilter("ignore")
            totals = []
            for k,v in summary.items():
                totals.append((k[0], k[1], v[0], v[0]+v[1], 100*v[0]/(v[0]+v[1]), np.nanmean(v[2])*1e6, np.nanmean(v[3])*1e6, np.nanmean(v[2]+v[3])*1e6))
            
        colsize = max([len(t[0]) + len(t[1]) for t in totals]) + 8
        totals.sort(key=lambda x: (x[4], x[3], x[0], x[1]), reverse=True)
        for tot in totals:
            pad = " " * (colsize - (len(tot[0]) + len(tot[1]) + 3))
            fields = list(tot)
            fields.insert(2, pad)
            fields = tuple(fields)
            try:
                print(u"%s → %s%s\t:\t%d/%d\t%0.2f%%\t%0.2f\t%0.2f\t%0.2f" % fields)
            except UnicodeEncodeError:
                print("%s - %s%s\t:\t%d/%d\t%0.2f%%\t%0.2f\t%0.2f\t%0.2f" % fields)

        print("\nTotal:  \t%d/%d\t%0.2f%%" % (tot_connected, tot_connected+tot_probed, 100*tot_connected/(tot_connected+tot_probed)))
        print("")

    def print_label_summary(self):
        print("-----------------------")
        print("       Labeling")
        print("-----------------------")

        n_qc_passed = 0
        n_dye_passed = 0
        n_qc_and_biocytin_passed = 0
        n_dye_and_biocytin_passed = 0

        for expt in self:
            for cell in expt.cells.values():
                biocytin = cell.labels.get('biocytin', None)
                if biocytin is None:
                    # ignore cells with no biocytin data
                    continue
                biocytin = (biocytin == '+')
                if cell.pass_qc:
                    n_qc_passed += 1
                    if biocytin:
                        n_qc_and_biocytin_passed += 1
                if cell.cre_type is not None:
                    n_dye_passed += 1
                    if biocytin:
                        n_dye_and_biocytin_passed += 1

        dye_biocytin_percent = 100*n_dye_and_biocytin_passed/n_dye_passed if n_dye_passed > 0 else 0
        print("%0.2f (%d/%d) of dye-filled cells had a biocytin fill" % (dye_biocytin_percent, n_dye_and_biocytin_passed, n_dye_passed))

        qc_biocytin_percent = 100*n_qc_and_biocytin_passed/n_qc_passed if n_qc_passed > 0 else 0
        print("%0.2f (%d/%d) of qc-passed cells had a biocytin fill" % (qc_biocytin_percent, n_qc_and_biocytin_passed, n_qc_passed))

        print("")


            
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--region', type=str)
    parser.add_argument('--start', type=arg_to_date)
    parser.add_argument('--stop', type=arg_to_date)
    parser.add_argument('--list-stims', action='store_true', default=False, dest='list_stims')
    parser.add_argument('files', nargs='+')
    args = parser.parse_args(sys.argv[1:])

    all_expts = ExperimentList()

    for f in args.files:
        all_expts.load(f)

    if len(all_expts) == 0:
        print("No experiments loaded; bailing out.")
        sys.exit(-1)

    expts = all_expts.select(start=args.start, stop=args.stop, region=args.region)

    if len(expts) == 0:
        print("No experiments selected; bailing out.")
        sys.exit(-1)

    expts.check()
    

    # Generate summary of experiments
    expts.print_expt_summary(args.list_stims)

    # Generate a summary of connectivity
    expts.print_connectivity_summary()

    # Print extra information about labeling
    expts.print_label_summary()


    p = pg.plot()
    expts.distance_plot('sim1', 'sim1', plot=p, color=(0, 0, 255))
    expts.distance_plot('tlx3', 'tlx3', plot=p, color=(200, 200, 0))
    expts.distance_plot('pvalb', 'pvalb', plot=p, color=(200, 0, 200))

    types = ['unknown', 'sim1', 'tlx3', 'pvalb', 'sst', 'vip']
    #types = ['sim1', 'unknown']
    expts.matrix(types, types)
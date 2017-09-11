# *-* coding: utf-8 *-*
from __future__ import print_function, division
import datetime
import glob
import json
import numpy as np
import os
import re

import pyqtgraph as pg
import pyqtgraph.configfile

from lims import specimen_info, specimen_images
from constants import ALL_CRE_TYPES, ALL_LABELS
from cell import Cell
from data import MultipatchExperiment


class Experiment(object):
    def __init__(self, entry):
        self.entry = entry
        self._connections = []
        self._region = None
        self._summary = None
        self._view = None
        self._site_info = None
        self._slice_info = None
        self._expt_info = None
        self._lims_record = None
        self._site_path = None
        self._probed = None
        self._sweep_summary = None
        self._mosaic_file = None
        self._nwb_file = None
        self._data = None
        self._stim_list = None

        try:
            self.source_id = self.id_from_entry(entry)
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
                raise Exception("Error parsing %s for experiment: %s\n%s" % (ch.lines[0], self, exc.args))

        # gather lists of all labels and cre types
        cre_types = set()
        labels = set()
        for cell in self.cells.values():
            if cell.cre_type not in cre_types:
                cre_types.add(cell.cre_type)
            labels |= set(cell.labels.keys()) & set(ALL_LABELS)
        self.cre_types = sorted(list(cre_types), key=lambda x: ALL_CRE_TYPES.index(x.split(',')[0]))
        self.labels = sorted(list(labels), key=lambda x: ALL_LABELS.index(x))

        # make sure all cells have information for all labels
        for cell in self.cells.values():
            for label in labels:
                assert label in cell.labels
            for cre in cre_types:
                for crepart in cre.split(','):
                    if crepart != 'unknown' and crepart not in cell.labels:
                        raise Exception('Cre type "%s" not in cell.labels: %s' % (crepart, cell.labels.keys()))

        # read cell positions from mosaic files
        try:
            self.load_cell_positions()
        except Exception as exc:
            print("Warning: Could not load cell positions for %s:\n    %s" % (self, exc.args[0]))

        # pull donor/specimen info from LIMS
        self.age
        # check for a single NWB file
        self.nwb_file

    @staticmethod
    def id_from_entry(entry):
        return (entry.file, entry.lines[0])

    @property
    def uid(self):
        """Return a unique ID string for this experiment.
        
        This returns the site timestamp formatted to 2 decimal places, which is
        very likely to be unique for any site.
        """
        return '%0.2f' % (self.site_info['__timestamp__'])

    @property
    def connections(self):
        """A list of connections reported for this experiment, excluding any that did not pass QC.
        
        Each item in the list is a tuple containing the pre- and postsynaptic cell IDs::
        
            [(pre_cell_id, post_cell_id), ...]
        """
        probed = self.connections_probed
        return [c for c in self._connections if c in probed]

    @property
    def sweep_summary(self):
        """A structure providing basic metadata on all sweeps collected in this
        experiment::
        
            [{dev_1: [stim_name, clamp_mode, holding_current, holding_potential], ...}, ...]
        """
        if self._sweep_summary is None:
            sweeps = []
            with self.data as nwb:
                for srec in nwb.contents:
                    sweep = {}
                    for dev in srec.devices:
                        rec = srec[dev]
                        sweep[dev] = rec.meta['stim_name'], rec.clamp_mode, rec.holding_current, rec.holding_potential
                    sweeps.append(sweep)
            self._sweep_summary = sweeps
        return self._sweep_summary

    def list_stims(self):
        """Open NWB file and return a list of stim set names.
        """
        if self._stim_list is None:
            stims = []
            for sweep in self.sweep_summary:
                for dev,info in sweep.items():
                    stim = info[0]
                    if stim not in stims:
                        stims.append(stim)

            # Shorten stim names
            stims = [self._short_stim_name(stim) for stim in stims]

            # sort by frequency
            def freq(stim):
                m = re.match('(.*)(\d+)Hz', stim)
                if m is None:
                    return (0, 0)
                else:
                    return (len(m.groups()[0]), int(m.groups()[1]))
            stims.sort(key=freq)

            self._stim_list = stims

        return self._stim_list

    @staticmethod
    def _short_stim_name(stim):
        if stim.startswith('PulseTrain_'):
            stim = stim[11:]
        elif stim.startswith('SPulseTrain_'):
            stim = 'S' + stim[12:]
        if stim.endswith('_DA_0'):
            stim = stim[:-5]
        if stim.endswith('H'):
            stim += 'z'
        return stim

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
            self._connections.append((int(m.groups()[0]), int(m.groups()[1])))

    def summary(self):
        """Return a structure summarizing (non)connectivity in the experiment.
        
        Looks like:
        
            {(pre_type, post_type): {
                'connected': n, 
                'unconnected': m, 
                'cdist': [...], 
                'udist': [...]
                }, 
            ...}
        """
        if self._summary is None:
            csum = {}
            for i, j in self.connections_probed:
                ci, cj = self.cells[i], self.cells[j]
                typ = (ci.cre_type, cj.cre_type)
                if typ not in csum:
                    csum[typ] = {'connected': 0, 'unconnected': 0, 'cdist':[], 'udist':[]}
                if (i, j) in self.connections:
                    csum[typ]['connected'] += 1
                    csum[typ]['cdist'].append(ci.distance(cj))
                else:
                    csum[typ]['unconnected'] += 1
                    csum[typ]['udist'].append(ci.distance(cj))
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
                        #print("Ignore probe (ind. cell type) %s:%d-%d" % (self.source_id, i, j))
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
        return sum([x['connected'] for x in self.summary().values()])

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
        if self._mosaic_file is None:
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
            self._mosaic_file = sitefile
        return self._mosaic_file

    @property
    def path(self):
        """Filesystem path to the root of this experiment.
        """
        if self._site_path is None:
            date, slice, site = self.source_id[1].split('-')
            root = os.path.dirname(self.source_id[0])
            if '_' not in date:
                date += '_000'
            paths = [
                os.path.join(root, date, "slice_%03d"%int(slice), "site_%03d"%int(site)),
                os.path.join(root, 'V1', date, "slice_%03d"%int(slice), "site_%03d"%int(site)),
                os.path.join(root, 'ALM', date, "slice_%03d"%int(slice), "site_%03d"%int(site)),
                # missing data, still in versioned backups
                os.path.join(root, '..', '..', '..', 'version_backups', 'data', 'Alex', 'V1', date, "slice_%03d" % int(slice), "site_%03d" % int(site)),
            ]
            for path in paths:
                if os.path.isdir(path):
                    self._site_path = path
                    break
            if self._site_path is None:
                raise Exception("Cannot find filesystem path for experiment %s. Attempted paths:\n%s" % (self, "\n".join(paths)))
        return self._site_path

    def __repr__(self):
        return "<Experiment %s (%s:%d)>" % (self.source_id[1], self.source_id[0], self.entry.lineno)

    @property
    def site_info(self):
        if self._site_info is None:
            index = os.path.join(self.path, '.index')
            if not os.path.isfile(index):
                raise TypeError("Cannot find slice index file (%s) for experiment %s" % (index, self))
            self._site_info = pg.configfile.readConfigFile(index)['.']
        return self._site_info

    @property
    def slice_info(self):
        if self._slice_info is None:
            index = os.path.join(os.path.split(self.path)[0], '.index')
            if not os.path.isfile(index):
                raise TypeError("Cannot find slice index file (%s) for experiment %s" % (index, self))
            self._slice_info = pg.configfile.readConfigFile(index)['.']
        return self._slice_info

    @property
    def expt_info(self):
        if self._expt_info is None:
            index = os.path.join(self.path, '..', '..', '.index')
            if not os.path.isfile(index):
                raise TypeError("Cannot find index file (%s) for experiment %s" % (index, self))
            self._expt_info = pg.configfile.readConfigFile(index)['.']
        return self._expt_info

    @property
    def nwb_file(self):
        if self._nwb_file is None:
            p = self.path
            files = glob.glob(os.path.join(p, '*.nwb'))
            if len(files) == 0:
                files = glob.glob(os.path.join(p, '*.NWB'))
            if len(files) == 0:
                raise Exception("No NWB file found for %s" % self)
            elif len(files) > 1:
                raise Exception("Multiple NWB files found for %s" % self)
            self._nwb_file = files[0]
        return self._nwb_file

    @property
    def data(self):
        """Data object from NWB file. 
        
        Contains all ephys recordings.
        """
        if self._data is None:
            if not os.path.isdir('cache'):
                os.mkdir('cache')
            cf = os.path.join('cache', self.nwb_file.replace('/', '_').replace(':', '_').replace('\\', '_'))
            if not os.path.isfile(cf) or os.stat(self.nwb_file).st_mtime > os.stat(cf).st_mtime:
                try:
                    import shutil
                    print("copying to cache:", cf)
                    shutil.copyfile(self.nwb_file, cf)
                except:
                    if os.path.isfile(cf):
                        os.remove(cf)
                    raise
            
            self._data = MultipatchExperiment(cf)
        return self._data

    def close_data(self):
        self.data.close()
        self._data = None

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
        """A dictionary of specimen information queried from LIMS.
        
        See multipatch_analysis.lims.section_info()
        """
        if self._lims_record is None:
            self._lims_record = specimen_info(self.specimen_id)
        return self._lims_record

    @property
    def biocytin_image_url(self):
        """A LIMS URL that points to the biocytin image for this specimen, or
        None if no image is found.
        """
        images = specimen_images(self.specimen_id)
        for img_id, treatment in images:
            if treatment == 'Biocytin':
                return "http://lims2/siv?sub_image=%d" % img_id

    @property
    def dapi_image_url(self):
        """A LIMS URL that points to the DAPI image for this specimen, or
        None if no image is found.
        """
        images = specimen_images(self.specimen_id)
        for img_id, treatment in images:
            if treatment == 'DAPI':
                return "http://lims2/siv?sub_image=%d" % img_id

    @property
    def multipatch_log(self):
        files = [p for p in os.listdir(self.path) if re.match(r'MultiPatch_\d+.log', p)]
        if len(files) == 0:
            raise TypeError("Could not find multipatch log file for %s" % self)
        if len(files) > 1:
            raise TypeError("Found multiple multipatch log files for %s" % self)
        return os.path.join(self.path, files[0])

    @property
    def surface_depth(self):
        try:
            mplog = self.multipatch_log
        except TypeError:
            return None
        lines = [l for l in open(mplog, 'rb').readlines() if 'surface_depth_changed' in l]
        if len(lines) == 0:
            return None
        line = lines[-1].rstrip(',\r\n')
        return json.loads(line)['surface_depth']

    @property
    def date(self):
        y,m,d = self.source_id[1].split('-')[0].split('_')[0].split('.')
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
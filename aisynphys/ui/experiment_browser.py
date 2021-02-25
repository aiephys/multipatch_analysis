from __future__ import print_function, division
from datetime import datetime
import pyqtgraph as pg

from aisynphys.database import default_db as db


class ExperimentBrowser(pg.TreeWidget):
    """TreeWidget showing a list of experiments with cells and pairs.
    """
    # TODO: add filtering options, cell info, context menu of actions, etc.
    
    def __init__(self):
        self.all_columns = ['date', 'timestamp', 'rig', 'organism', 'project', 'region', 'genotype', 'acsf']
        self.visible_columns = self.all_columns[:]
        
        pg.TreeWidget.__init__(self)
        
        self.setColumnCount(len(self.all_columns))
        self.setHeaderLabels(self.all_columns)
        self.setDragDropMode(self.NoDragDrop)
        self._last_expanded = None
        
    def populate(self, experiments=None, all_pairs=False, synapses=False):
        """Populate the browser with a list of experiments.
        
        Parameters
        ----------
        experiments : list | None
            A list of Experiment instances. If None, then automatically query experiments from the default database.
        all_pairs : bool
            If False, then pairs with no qc-passed pulse responses are excluded
        synapses : bool
            If True, then only synaptically connected pairs are shown
        """
        with pg.BusyCursor():
            # if all_pairs is set to True, all pairs from an experiment will be included regardless of whether they have data
            self.items_by_pair_id = {}
            
            self.session = db.session()
            
            if experiments is None:
                # preload all cells,pairs so they are not queried individually later on
                q = self.session.query(db.Pair, db.Experiment, db.Cell, db.Slice)
                q = q.join(db.Experiment, db.Pair.experiment_id==db.Experiment.id)
                q = q.join(db.Cell, db.Cell.id==db.Pair.pre_cell_id)
                q = q.join(db.Slice)
                if synapses:
                    q = q.filter(db.Pair.has_synapse==True)
                    
                recs = q.all()
                experiments = list(set([rec.Experiment for rec in recs]))
            
            experiments.sort(key=lambda e: e.acq_timestamp)
            for expt in experiments:
                date = expt.acq_timestamp
                date_str = datetime.fromtimestamp(date).strftime('%Y-%m-%d')
                slice = expt.slice
                expt_item = pg.TreeWidgetItem(map(str, [date_str, '%0.3f'%expt.acq_timestamp, expt.rig_name, slice.species, expt.project_name, expt.target_region, slice.genotype, expt.acsf]))
                expt_item.expt = expt
                self.addTopLevelItem(expt_item)

                for pair in expt.pair_list:
                    if all_pairs is False and pair.n_ex_test_spikes == 0 and pair.n_in_test_spikes == 0:
                        continue
                    if synapses and not pair.has_synapse:
                        continue
                    cells = '%s => %s' % (pair.pre_cell.ext_id, pair.post_cell.ext_id)
                    conn = {True:"syn", False:"-", None:"?"}[pair.has_synapse]
                    gap = {True:"gap", False:"-", None:"?"}[pair.has_electrical]
                    pre_type = pair.pre_cell.cre_type
                    if pre_type == 'unknown':
                        dendrite_type = pair.pre_cell.morphology.dendrite_type if pair.pre_cell.morphology is not None else None
                        pre_type = dendrite_type if dendrite_type in ['spiny', 'aspiny', 'sparsely spiny'] else pre_type
                    post_type = pair.post_cell.cre_type
                    if post_type == 'unknown':
                        dendrite_type = pair.post_cell.morphology.dendrite_type if pair.pre_cell.morphology is not None else None
                        post_type = dendrite_type if dendrite_type in ['spiny', 'aspiny', 'sparsely spiny'] else post_type
                    types = 'L%s %s => L%s %s' % (pair.pre_cell.target_layer or "?", pre_type, pair.post_cell.target_layer or "?", post_type)
                    pair_item = pg.TreeWidgetItem([cells, conn+' ; '+gap, types])
                    expt_item.addChild(pair_item)
                    pair_item.pair = pair
                    pair_item.expt = expt
                    self.items_by_pair_id[pair.id] = pair_item
                    # also allow select by ext id
                    self.items_by_pair_id[(expt.acq_timestamp, pair.pre_cell.ext_id, pair.post_cell.ext_id)] = pair_item
                    
            self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())
                
    def select_pair(self, pair_id):
        """Select a specific pair from the list
        """
        if self._last_expanded is not None:
            self._last_expanded.setExpanded(False)
        item = self.items_by_pair_id[pair_id]
        self.clearSelection()
        item.setSelected(True)
        parent = item.parent()
        if not parent.isExpanded():
            parent.setExpanded(True)
            self._last_expanded = parent
        else:
            self._last_expanded = None
        self.scrollToItem(item)

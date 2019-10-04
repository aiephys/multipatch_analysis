import sys, argparse
import pyqtgraph as pg
from aisynphys.database import default_db as db
from aisynphys.matrix_analyzer import MatrixAnalyzer
from collections import OrderedDict

if __name__ == '__main__':

    app = pg.mkQApp()
    if sys.flags.interactive == 1:
        pg.dbg()
    # pg.setConfigOption('background', 'w')
    # pg.setConfigOption('foreground', 'k')

    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', type=str, default='external')
    args = parser.parse_args(sys.argv[1:])
    analyzer_mode = args.mode

    session = db.session()
    
    # Define cell classes
    cell_class_groups = OrderedDict([
        ('All Transgenic Classes', [
            # {'cre_type': 'unknown', 'target_layer': '2/3','cortical_layer': '2/3'},
            {'dendrite_type': 'spiny', 'target_layer': '2/3', 'cortical_layer': '2/3','display_names': ('L2/3', 'Pyr\nspiny')},
            {'cre_type': 'pvalb', 'target_layer': '2/3', 'cortical_layer': '2/3', 'display_names': ('L2/3', 'Pv')},
            {'cre_type': 'sst', 'target_layer': '2/3', 'cortical_layer': '2/3', 'display_names': ('L2/3', 'Sst')},
            {'cre_type': 'vip', 'target_layer': '2/3', 'cortical_layer': '2/3', 'display_names': ('L2/3', 'Vip')},
           # {'cre_type': 'rorb', 'target_layer': '4'},
            {'cre_type': 'nr5a1', 'target_layer': '4', 'cortical_layer': '4', 'display_names': ('L4', 'Pyr\n nr5a1')},
            {'cre_type': 'pvalb', 'target_layer': '4', 'cortical_layer': '4', 'display_names': ('L4', 'Pv')},
            {'cre_type': 'sst', 'target_layer': '4', 'cortical_layer': '4', 'display_names': ('L4', 'Sst')},
            {'cre_type': 'vip', 'target_layer': '4', 'cortical_layer': '4', 'display_names': ('L4', 'Vip')},
            {'cre_type': ('sim1', 'fam84b'), 'target_layer': '5', 'cortical_layer': '5', 'display_names': ('L5', 'Pyr ET\nsim1, fam84b')},
            {'cre_type': 'tlx3', 'target_layer': '5', 'display_names': ('L5', 'Pyr IT\ntlx3'), 'cortical_layer': '5'},
            {'cre_type': 'pvalb', 'target_layer': '5', 'cortical_layer': '5', 'display_names': ('L5', 'Pv')},
            {'cre_type': 'sst', 'target_layer': '5', 'cortical_layer': '5', 'display_names': ('L5', 'Sst')},
            {'cre_type': 'vip', 'target_layer': '5', 'cortical_layer': '5', 'display_names': ('L5', 'Vip')},
            {'cre_type': 'ntsr1', 'target_layer': '6', 'cortical_layer': ('6a', '6b'), 'display_names': ('L6', 'Pyr\nntsr1')},
            {'cre_type': 'pvalb', 'target_layer': '6', 'cortical_layer': ('6a', '6b'), 'display_names': ('L6', 'Pv')},
            {'cre_type': 'sst', 'target_layer': '6', 'cortical_layer': ('6a', '6b'), 'display_names': ('L6', 'Sst')},
            {'cre_type': 'vip', 'target_layer': '6', 'cortical_layer': ('6a', '6b'), 'display_names': ('L6', 'Vip')},
        ]),

        ('Mouse Layer 2/3', [
            # {'cre_type': 'unknown', 'target_layer': '2/3', 'cortical_layer': '2/3'},
            #{'pyramidal': True, 'target_layer': '2/3'},
            {'dendrite_type': 'spiny', 'target_layer': '2/3', 'cortical_layer': '2/3'},
            {'cre_type': 'pvalb', 'target_layer': '2/3', 'cortical_layer': '2/3'},
            {'cre_type': 'sst', 'target_layer': '2/3', 'cortical_layer': '2/3'},
            {'cre_type': 'vip', 'target_layer': '2/3', 'cortical_layer': '2/3'},
        ]),
        
        ('Mouse Layer 4', [
            {'cre_type': 'nr5a1', 'target_layer': '4', 'cortical_layer': '4'},
            {'cre_type': 'pvalb', 'target_layer': '4', 'cortical_layer': '4'},
            {'cre_type': 'sst', 'target_layer': '4', 'cortical_layer': '4'},
            {'cre_type': 'vip', 'target_layer': '4', 'cortical_layer': '4'},
        ]),

        ('Mouse Layer 5', [
            {'cre_type': ('sim1', 'fam84b'), 'target_layer': '5', 'display_names': ('L5', 'PT\nsim1, fam84b'), 'cortical_layer': '5'},
            {'cre_type': 'tlx3', 'target_layer': '5', 'display_names': ('L5', 'IT\ntlx3'), 'cortical_layer': '5'},
            {'cre_type': 'pvalb', 'target_layer': '5', 'cortical_layer': '5'},
            {'cre_type': 'sst', 'target_layer': '5', 'cortical_layer': '5'},
            {'cre_type': 'vip', 'target_layer': '5', 'cortical_layer': '5'},
        ]),

        ('Mouse Layer 6', [
            {'cre_type': 'ntsr1', 'target_layer': '6', 'cortical_layer': ('6a', '6b'), 'display_names': ('L6', 'ntsr1')},
            {'cre_type': 'pvalb', 'target_layer': '6', 'cortical_layer': ('6a', '6b'), 'display_names': ('L6', 'pvalb')},
            {'cre_type': 'sst', 'target_layer': '6', 'cortical_layer': ('6a', '6b'), 'display_names': ('L6', 'sst')},
            {'cre_type': 'vip', 'target_layer': '6', 'cortical_layer': ('6a', '6b'), 'display_names': ('L6', 'vip')},
        ]),

        ('Inhibitory Transgenic Classes',[
            {'cre_type': 'pvalb', 'display_names': ('', 'Pv')},
            {'cre_type': 'sst', 'display_names': ('', 'Sst')},
            {'cre_type': 'vip', 'display_names': ('', 'Vip')},
        ]),
 
        ('Excitatory Transgenic Classes', [
            {'dendrite_type': 'spiny', 'target_layer': '2/3', 'cortical_layer': '2/3', 'display_names': ('L2/3', 'Pyr\nspiny')},
            # {'cre_type': 'unknown', 'target_layer': '2/3'},
            {'cre_type': 'nr5a1', 'target_layer': '4', 'cortical_layer': '4', 'display_names': ('L4', 'Pyr\nnr5a1')},
            {'cre_type': ('sim1', 'fam84b'), 'target_layer': '5', 'display_names': ('L5', 'Pyr ET\nsim1, fam84b'), 'cortical_layer': '5'},
            {'cre_type': 'tlx3', 'target_layer': '5', 'display_names': ('L5', 'Pyr IT\ntlx3'), 'cortical_layer': '5'},
            {'cre_type': 'ntsr1', 'target_layer': '6', 'cortical_layer': ('6a', '6b'), 'display_names': ('L6', 'Pyr\nntsr1')}
        ]),

        ('Mouse E-I Cre-types', [
            {'cre_type': ('unknown', 'nr5a1', 'tlx3', 'sim1', 'ntsr1'), 'display_names': ('', 'Excitatory\nunknown, nr5a1,\ntlx3, sim1, ntsr1')},
            {'cre_type': ('pvalb', 'sst', 'vip'), 'display_names': ('', 'Inhibitory\npvalb, sst, vip')},
        ]),

        ('Inhibitory Transgenic Classes by layer',[
            # {'pyramidal': True, 'target_layer': '2/3'},
            # {'cre_type': 'unknown', 'target_layer': '2/3'},
            {'cre_type': ('pvalb', 'sst', 'vip'), 'target_layer': '2/3', 'display_names': ('L2/3', 'Inhibitory\npv, sst, vip'), 'cortical_layer': '2/3'},
            # {'cre_type': 'nr5a1', 'target_layer': '4'},
            {'cre_type': ('pvalb', 'sst', 'vip'), 'target_layer': '4', 'display_names': ('L4', 'Inhibitory\npv, sst, vip'), 'cortical_layer': '4'},
            # {'cre_type': 'sim1', 'target_layer': '5'},
            # {'cre_type': 'tlx3', 'target_layer': '5'},
            {'cre_type': ('pvalb', 'sst', 'vip'), 'target_layer': '5', 'display_names': ('L5', 'Inhibitory\npv, sst, vip'), 'cortical_layer': '5'},
            # {'cre_type': 'ntsr1', 'target_layer': '6'},
            {'cre_type': ('pvalb', 'sst', 'vip'), 'target_layer': '6', 'display_names': ('L6', 'Inhibitory\npv, sst, vip'), 'cortical_layer': ('6a', '6b')},     
        ]),

        ('Pyramidal / Nonpyramidal by layer', [
            {'pyramidal': True, 'target_layer': '2'},
            {'pyramidal': False, 'target_layer': '2'},
            {'pyramidal': True, 'target_layer': '3'},
            {'pyramidal': False, 'target_layer': '3'},
            {'pyramidal': True, 'target_layer': '4'},
            {'pyramidal': False, 'target_layer': '4'},
            {'pyramidal': True, 'target_layer': '5'},
            {'pyramidal': False, 'target_layer': '5'},
            {'pyramidal': True, 'target_layer': '6'},
            {'pyramidal': False, 'target_layer': '6'},
        ]),

        ('Pyramidal Cells', [
            {'dendrite_type': 'spiny', 'target_layer': '2', 'cortical_layer': '2', 'display_names': ('L2', 'Pyr\nspiny')},
            {'dendrite_type': 'spiny', 'target_layer': '2/3', 'cortical_layer': '2/3', 'display_names': ('L2/3', 'Pyr\nspiny')}, 
            {'dendrite_type': 'spiny', 'target_layer': '3', 'cortical_layer': '3', 'display_names': ('L3', 'Pyr\nspiny')},
            {'dendrite_type': 'spiny', 'target_layer': '4', 'cortical_layer': '4', 'display_names': ('L4', 'Pyr\nspiny')},
            {'dendrite_type': 'spiny', 'target_layer': '5', 'cortical_layer': '5', 'display_names': ('L5', 'Pyr\nspiny')},
            {'dendrite_type': 'spiny', 'target_layer': '6','cortical_layer': ('6','6a', '6b'), 'display_names': ('L6', 'Pyr\nspiny')},
        ]),

        ('Non-Pyramidal Cells', [
            {'dendrite_type': 'aspiny', 'target_layer': '2', 'cortical_layer': '2', 'display_names': ('L2', 'Non-Pyr\naspiny')},
            {'dendrite_type': 'aspiny', 'target_layer': '2/3', 'cortical_layer': '2/3', 'display_names': ('L2/3', 'Non-Pyr\naspiny')}, 
            {'dendrite_type': 'aspiny', 'target_layer': '3', 'cortical_layer': '3', 'display_names': ('L3', 'Non-Pyr\naspiny')},
            {'dendrite_type': 'aspiny', 'target_layer': '4', 'cortical_layer': '4', 'display_names': ('L4', 'Non-Pyr\naspiny')},
            {'dendrite_type': 'aspiny', 'target_layer': '5', 'cortical_layer': '5', 'display_names': ('L5', 'Non-Pyr\naspiny')},
            {'dendrite_type': 'aspiny', 'target_layer': '6','cortical_layer': ('6', '6a', '6b'), 'display_names': ('L6', 'Non-Pyr\naspiny')},
        ]),

        ('All Cells', [
            {'target_layer': '2', 'cortical_layer': '2'},
            {'target_layer': '2/3', 'cortical_layer': '2/3'},
            {'target_layer': '3', 'cortical_layer': '3'},
            {'target_layer': '4', 'cortical_layer': '4'},
            {'target_layer': '5', 'cortical_layer': '5'},
            {'target_layer': '6', 'cortical_layer': ('6', '6a', '6b'), 'display_names': ('', 'L6')},
        ]),

        ('2P-Opto cre types', [
            {'cre_type':'ntsr1'},
            #{'cre_type':'unknown'},
            {'cre_type':'sst'},
            {'cre_type':'tlx3'},
            {'cre_type':'rorb'},
            {'cre_type':'scnn1a'}])
    ])

    if analyzer_mode == 'external':
        groups = ['All Transgenic Classes','Excitatory Transgenic Classes', 'Inhibitory Transgenic Classes', 'Inhibitory Transgenic Classes by layer', 'All Cells', 'Pyramidal Cells', 'Non-Pyramidal Cells']
        cell_class_groups = {g:cell_class_groups[g] for g in groups}

    maz = MatrixAnalyzer(session=session, cell_class_groups=cell_class_groups, default_preset='None', preset_file='matrix_analyzer_presets.json', analyzer_mode=analyzer_mode)

    if sys.flags.interactive == 0:
        app.exec_()
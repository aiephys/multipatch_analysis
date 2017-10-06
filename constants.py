
ALL_LABELS = ['biocytin', 'af488', 'cascade_blue']
INHIBITORY_CRE_TYPES = ['sst','pvalb', 'vip', 'ndnf', 'chat', 'htr3a', 'nos1', 'chrna2']
EXCITATORY_CRE_TYPES = ['tlx3', 'sim1', 'rorb', 'ntsr1', 'rbp4', 'ctgf', 'glt25d2', 'L23pyr']
ALL_CRE_TYPES = INHIBITORY_CRE_TYPES + EXCITATORY_CRE_TYPES + ['unknown']


REPORTER_LINES = {
    'Ai2':                        'EYFP',
    'Ai3':                        'EYFP',
    'Ai6':                        'ZsGreen',
    'Ai9':                        'tdTomato',
    'Ai14':                       'tdTomato',
    'Ai65(RCFL-tdT)':             'tdTomato',
    'Ai65F':                      'tdTomato',
    'Ai66(RCRL-tdT)':             'tdTomato',
    'Ai57(RCFL-Jaws)':            'Jaws-GFP-ER2',
    'Ai72(RCL-VSFPB)':            'VSFP-Butterfly 1.2',
    'Ai87(RCL-iGluSnFR)':         'iGluSnFR',
    'Ai95(RCL-GCaMP6f)':          'GCaMP6f',
    'Ai96(RCL-GCaMP6s)':          'GCaMP6s',
    'Ai62(TITL-tdT)':             'tdTomato',
    'Ai82(TITL-GFP)':             'EGFP',
    'Ai79(TITL-Jaws)':            'Jaws-GFP-ER2',
    'Ai93(TITL-GCaMP6f)':         'GCaMP6f',
    'Ai94(TITL-GCaMP6s)':         'GCaMP6s',
    'Ai92(TITL-YCX2.60)':         'YCX2.60',
    'Ai78(TITL-VSFPB)':           'VSFP-Butterfly 1.2',
    'Ai85(TITL-iGluSnFR)':        'iGluSnFR',
    'Ai139(TIT2L-GFP-ICL-TPT)-D': 'EGFP + tdTomato',
}
   

FLUOROPHORES = {
    'tdTomato': 'red',
    'GFP': 'green',
    'EGFP': 'green',
    'AF488': 'green',
    'Cascade Blue': 'blue',
    'EYFP': 'yellow',
    'ZsGreen': 'green',
}



ACSF_RECIPES = ["2mM Ca & Mg", "1.3mM Ca & 1mM Mg"]
INTERNAL_RECIPES = ["Standard K-Gluc"]
INTERNAL_DYES = ['Cascade Blue', 'AF488']


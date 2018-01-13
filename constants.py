ALL_LABELS = ['biocytin', 'af488', 'cascade_blue']
HUMAN_LABELS = ['human_L2', 'human_L3', 'human_L4', 'human_L5', 'human_L6']
INHIBITORY_CRE_TYPES = ['sst','pvalb', 'vip', 'ndnf', 'chat', 'htr3a', 'nos1', 'chrna2']
EXCITATORY_CRE_TYPES = ['tlx3', 'sim1', 'rorb', 'ntsr1', 'rbp4', 'ctgf', 'glt25d2', 'L23pyr', 'slc17a8'] + HUMAN_LABELS
ALL_CRE_TYPES = INHIBITORY_CRE_TYPES + EXCITATORY_CRE_TYPES + ['unknown']


# see: 
# http://help.brain-map.org/download/attachments/2818171/Connectivity_Resources.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4365051/pdf/nihms665520.pdf

# RCL  : Rosa_CAG-LSL      (cre dependent)
# RCFL : Rosa_CAG-FSF-LSL  (cre+flp dependent)
# RCRL : Rosa_CAG-RSR-LSL  (cre+dre dependent)
# TIT : TIGRE-insulators-TRE_promoter   (tTa|rtTA dependent)
# TITL : TIGRE-insulators-TRE_promoter-LSL   (cre+(tTA|rtTA) dependent)
# TIT2L : TIGRE-insulators-TRE_promoter-LSL   (cre+(tTA|rtTA) dependent)
# IRES :
# T2A :

# LSL : lox-stop-lox (cre dependent)
# FSF : frt-stop-frt (flp dependent)
# RSR : rox-stop-rox (dre dependent, but also slightly activated by cre)

# tTA : tet-off  (activates only in absence of doxycycline)
# rtTA: tet-on   (activates only in presence of doxycycline)
# TRE : tetracycline response element - series of TetO sequences that
#       require tTA or rtTA to increase downstream expression
#         - dox binds tTa, which prevents further binding to TetO
#         - dox binds tTa, which allows further binding to TetO


REPORTER_LINES = {
    'Ai2':                        ('',           'EYFP'),
    'Ai3':                        ('',           'EYFP'),
    'Ai6':                        ('',           'ZsGreen'),
    'Ai9':                        ('',           'tdTomato'),
    'Ai14':                       ('',           'tdTomato'),
    'Ai14(RCL-tdT)':              ('?',          'tdTomato?'),
    'Ai27':                       ('',           'hChR2(H134R)-tdTomato'),
    'Ai31':                       ('',           'Syp-Emerald'),
    'Ai32':                       ('',           'ChR2(H134R)-EYFP'),
    'Ai34':                       ('',           'Syp-tdTomato'),
    'Ai35':                       ('',           'Arch-EGFP-ER2'),
    'Ai57(RCFL-Jaws)':            ('cre&flp',    'Jaws-GFP-ER2'),
    'Ai62(TITL-tdT)':             ('cre&tTA',    'tdTomato'),
    'Ai63(TIT-tdT)':              ('tTA',        'tdTomato'),
    'Ai65(RCFL-tdT)':             ('cre&flp',    'tdTomato'),
    'Ai65F':                      ('flp',        'tdTomato'),
    'Ai66(RCRL-tdT)':             ('cre&dre',    'tdTomato'),
    'Ai72(RCL-VSFPB)':            ('cre',        'VSFP-Butterfly 1.2'),
    'Ai87(RCL-iGluSnFR)':         ('cre',        'iGluSnFR'),
    'Ai95(RCL-GCaMP6f)':          ('cre',        'GCaMP6f'),
    'Ai96(RCL-GCaMP6s)':          ('cre',        'GCaMP6s'),
    'Ai62(TITL-tdT)':             ('cre&tTA',    'tdTomato'),
    'Ai82(TITL-GFP)':             ('cre&tTA',    'EGFP'),
    'Ai79(TITL-Jaws)':            ('cre&tTA',    'Jaws-GFP-ER2'),
    'Ai93(TITL-GCaMP6f)':         ('cre&tTA',    'GCaMP6f'),
    'Ai94(TITL-GCaMP6s)':         ('cre&tTA',    'GCaMP6s'),
    'Ai92(TITL-YCX2.60)':         ('cre&tTA',    'YCX2.60'),
    'Ai78(TITL-VSFPB)':           ('cre&tTA',    'VSFP-Butterfly 1.2'),
    'Ai85(TITL-iGluSnFR)':        ('cre&tTA',    'iGluSnFR'),
    'Ai139(TIT2L-GFP-ICL-TPT)-D': ('cre',        'EGFP+tdTomato'),
    'Ai140(TIT2L-GFP-ICL-tTA2)':  ('cre',        'EGFP'),
}
   
   
DRIVER_LINES = {
    'Rorb-T2A-tTA2':              ('cre,tTA',    'rorb'),
    'Tlx3-Cre_PL56':              ('cre',        'tlx3'),
    'Sim1-Cre_KJ18':              ('cre',        'sim1'),
    'Ntsr1-Cre_GN220':            ('cre',        'ntsr1'),
    'Chat-IRES-Cre-neo':          ('cre',        'chat'),
    'Rbp4-Cre_KL100':             ('cre',        'rbp4'),
    'Pvalb-IRES-Cre':             ('cre',        'pvalb'),
    'Sst-IRES-Cre':               ('cre',        'sst'),
    'Vip-IRES-Cre':               ('cre',        'vip'),
    'Pvalb-T2A-FlpO':             ('flp',        'pvalb'),
    'Sst-IRES-FlpO':              ('flp',        'sst'),
    'Vip-IRES-FlpO':              ('flp',        'vip'),
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
INHIBITORY_CRE_TYPES = ['sst', 'pvalb', 'vip', 'ndnf', 'chat', 'htr3a', 'nos1', 'chrna2']
EXCITATORY_CRE_TYPES = ['tlx3', 'sim1', 'rorb', 'ntsr1', 'rbp4', 'ctgf', 'glt25d2', 'slc17a8', 'cux2', 'nr5a1']
ALL_CRE_TYPES = INHIBITORY_CRE_TYPES + EXCITATORY_CRE_TYPES + ['unknown']

DRIVER_LAYERS = {
    'cux2':    ['2', '2/3', '3', '4'],
    'rorb':    ['4'],
    'nr5a1':   ['4'],
    'tlx3':    ['5'],
    'sim1':    ['5'],
    'rbp4':    ['5'],
    'slc17a8': ['5', '6'],
    'ntsr1':   ['6'],
    'ctgf':    ['6'],
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
# note: the lower-case af488 and cascade_blue are for backward compatibility; these
# may be removed at some point
ALL_LABELS = ['biocytin', 'af488', 'cascade_blue'] + FLUOROPHORES.keys()


LAYERS = ['1', '2', '2/3', '3', '4', '5', '5a', '5b', '6']

ACSF_RECIPES = ["2mM Ca & Mg", "1.3mM Ca & 1mM Mg"]
INTERNAL_RECIPES = ["Standard K-Gluc"]
INTERNAL_DYES = ['Cascade Blue', 'AF488']


# ----------------- Genotype Handling ----------------------------
# The purpose of the code below is to be able to ask, given a mouse genotype string, 
# what cell types are marked, and what reporter color are they marked with.

# There are two approaches:
# 1. Write down a list of all possible genotype strings and the associated 
#    driver / reporter information. This is simple and transparent, but requires
#    constant maintenance as we add more genotypes.
# 2. Write down the complete list of drivers and reporters, and use some simple
#    genetic modeling to predict what phenotype will be created by a particular
#    genotype.
# Method (2) is mostly implemented below, but encounters some complexity with the "genetic
# modeling" aspect. This is surely resolvable, but for now we just use (1) for its
# simplicity.

# -------------------------------------------------------------------------------
# Here is the simple method (1)

# Each genotype string maps to a list of (driver,[reporters]) pairs
from collections import OrderedDict
GENOTYPES = OrderedDict([
    ("Tlx3-Cre_PL56/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'pvalb': ['tdTomato'], 'tlx3': ['EGFP']}),
    ("Vip-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'pvalb': ['tdTomato'], 'vip': ['EGFP']}),
    ("Ntsr1-Cre_GN220/wt;Pvalb-T2A-FlpO/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'ntsr1': ['EGFP']}),
    ("Tlx3-Cre_PL56/wt;Vip-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'vip': ['tdTomato'], 'tlx3': ['EGFP']}),
    ("Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'pvalb': ['tdTomato'], 'sst': ['EGFP']}),
    ("Tlx3-Cre_PL56/wt;Vip-IRES-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'vip': ['tdTomato'], 'tlx3': ['EGFP']}),
    ("Tlx3-Cre_PL56/wt;Pvalb-T2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'pvalb': ['tdTomato'], 'tlx3': ['EGFP']}),
    ("Ntsr1-Cre_GN220/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'pvalb': ['tdTomato'], 'ntsr1': ['EGFP']}),
    ("Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'pvalb': ['tdTomato'], 'sst': ['EGFP']}),
    ("Tlx3-Cre_PL56/wt;Vip-IRES-FlpO/wt;PhiC31-neo/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'vip': ['tdTomato'], 'tlx3': ['EGFP']}),
    ("Sim1-Cre_KJ18/wt;Sst-IRES-FlpO/wt;Ai65F/Ai65F;Ai139(TIT2L-GFP-ICL-TPT)/wt", {'sst': ['tdTomato'], 'sim1': ['EGFP', 'tdTomato']}),
    ("Tlx3-Cre_PL56/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'sst': ['tdTomato'], 'tlx3': ['EGFP']}),
    ("Vip-IRES-Cre/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'sst': ['tdTomato'], 'vip': ['EGFP']}),
    ("Sim1-Cre_KJ18/wt;Pvalb-T2A-FlpO/wt;Snap25-LSL-F2A-GFP/wt;Ai65F/wt", {'pvalb': ['tdTomato'], 'sim1': ['EGFP']}),
    ("Vip-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'pvalb': ['tdTomato'], 'vip': ['EGFP']}),
    ("Sim1-Cre_KJ18/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai139(TIT2L-GFP-ICL-TPT)/wt", {'sst': ['tdTomato'], 'sim1': ['EGFP', 'tdTomato']}),
    ("Tlx3-Cre_PL56/wt;Sst-IRES-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'sst': ['tdTomato'], 'tlx3': ['EGFP']}),
    ("Ntsr1-Cre_GN220/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt;Ai65F/wt", {'ntsr1': ['EGFP']}),
    ("Slc17a8-IRES2-Cre/wt;Ai14(RCL-tdT)/wt", {'slc17a8': ['tdTomato']}),
    ("Ntsr1-Cre_GN220/wt;Vip-IRES-FlpO/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt;Ai65F/wt", {'vip': ['tdTomato'], 'ntsr1': ['EGFP']}),
    ("Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/Ai65F", {'pvalb': ['tdTomato']}),
    ("Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt", {'pvalb': ['tdTomato']}),
    ("Ntsr1-Cre_GN220/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'sst': ['tdTomato'], 'ntsr1': ['EGFP']}),
    ("Sim1-Cre_KJ18/wt;Ai14(RCL-tdT)/wt", {'sim1': ['tdTomato']}),
    ("Pvalb-IRES-Cre/wt;Ai14(RCL-tdT)/wt", {'pvalb': ['tdTomato']}),
    ("Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'sst': ['tdTomato']}),
    ("Nr5a1-Cre/wt;Ai14(RCL-tdT)/wt", {'nr5a1': ['tdTomato']}),
    ("Tlx3-Cre_PL56/wt;Vip-IRES-FlpO/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'tlx3': ['EGFP']}),
    ("Ntsr1-Cre_GN220/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'ntsr1': ['EGFP']}),
    ("Sim1-Cre_KJ18/wt;Sst-IRES-FlpO/wt;Ai139(TIT2L-GFP-ICL-TPT)/wt", {'sim1': ['EGFP', 'tdTomato']}),
    ("Tlx3-Cre_PL56/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt", {'tlx3': ['EGFP']}),
    ("Pvalb-IRES-Cre/wt;Rorb-T2A-tTA2/wt;Ai63(TIT-tdT)/Ai140(TIT2L-GFP-ICL-tTA2)", {'pvalb': ['EGFP'], 'rorb': ['tdTomato']}),
])    




# -------------------------------------------------------------------------------
# Here is the more general/complex method (2) that attempts to automatically parse
# genotype strings (see also genotypes.py)


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

# TRE : tetracycline response element - series of TetO sequences that
#       require tTA or rtTA to increase downstream expression
#         - dox binds tTa, which prevents further binding to TetO
#         - dox binds rtTa, which allows further binding to TetO
#       this system is also used for amplification.
# tTA : tet-off  (activates TRE only in absence of doxycycline)
# rtTA: tet-on   (activates TRE only in presence of doxycycline)
# tTA2 : same behavior as tTA


EXPRESSION_FACTORS = ['cre', 'flp', 'dre', 'tTA']

REPORTER_LINES = {                # dependencies # products
    'Ai2':                        ('',           ['EYFP']),
    'Ai3':                        ('',           ['EYFP']),
    'Ai6':                        ('',           ['ZsGreen']),
    'Ai9':                        ('',           ['tdTomato']),
    'Ai14':                       ('',           ['tdTomato']),
    'Ai14(RCL-tdT)':              ('cre',        ['tdTomato']),
    'Ai27':                       ('',           ['hChR2(H134R)','tdTomato']),
    'Ai31':                       ('',           ['Syp','Emerald']),
    'Ai32':                       ('',           ['ChR2(H134R)','EYFP']),
    'Ai34':                       ('',           ['Syp','tdTomato']),
    'Ai35':                       ('',           ['Arch','EGFP','ER2']),
    'Ai57(RCFL-Jaws)':            ('cre&flp',    ['Jaws','GFP','ER2']),
    'Ai62(TITL-tdT)':             ('cre&tTA',    ['tdTomato']),
    'Ai63(TIT-tdT)':              ('tTA',        ['tdTomato']),
    'Ai65(RCFL-tdT)':             ('cre&flp',    ['tdTomato']),
    'Ai65F':                      ('flp',        ['tdTomato']),
    'Ai66(RCRL-tdT)':             ('cre&dre',    ['tdTomato']),
    'Ai72(RCL-VSFPB)':            ('cre',        ['VSFP','Butterfly 1.2']),
    'Ai87(RCL-iGluSnFR)':         ('cre',        ['iGluSnFR']),
    'Ai95(RCL-GCaMP6f)':          ('cre',        ['GCaMP6f']),
    'Ai96(RCL-GCaMP6s)':          ('cre',        ['GCaMP6s']),
    'Ai62(TITL-tdT)':             ('cre&tTA',    ['tdTomato']),
    'Ai82(TITL-GFP)':             ('cre&tTA',    ['EGFP']),
    'Ai79(TITL-Jaws)':            ('cre&tTA',    ['Jaws','GFP','ER2']),
    'Ai93(TITL-GCaMP6f)':         ('cre&tTA',    ['GCaMP6f']),
    'Ai94(TITL-GCaMP6s)':         ('cre&tTA',    ['GCaMP6s']),
    'Ai92(TITL-YCX2.60)':         ('cre&tTA',    ['YCX2.60']),
    'Ai78(TITL-VSFPB)':           ('cre&tTA',    ['VSFP','Butterfly 1.2']),
    'Ai85(TITL-iGluSnFR)':        ('cre&tTA',    ['iGluSnFR']),
    'Ai139(TIT2L-GFP-ICL-TPT)-D': ('cre',        ['EGFP','tdTomato']),
    'Ai139(TIT2L-GFP-ICL-TPT)':   ('cre',        ['EGFP','tdTomato']),
    'Ai140(TIT2L-GFP-ICL-tTA2)':  ('cre',        ['EGFP', 'tTA']),
    'Snap25-LSL-F2A-GFP':         ('cre',        ['EGFP']),
}
   
DRIVER_LINES = {                  # dependencies   # products
    'Nr5a1-Cre':                  ('nr5a1',        ['cre']),
    'Rorb-T2A-tTA2':              ('rorb',         ['tTA']),
    'Tlx3-Cre_PL56':              ('tlx3',         ['cre']),
    'Sim1-Cre_KJ18':              ('sim1',         ['cre']),
    'Ntsr1-Cre_GN220':            ('ntsr1',        ['cre']),
    'Chat-IRES-Cre-neo':          ('chat',         ['cre']),
    'Rbp4-Cre_KL100':             ('rbp4',         ['cre']),
    'Pvalb-IRES-Cre':             ('pvalb',        ['cre']),
    'Sst-IRES-Cre':               ('sst',          ['cre']),
    'Vip-IRES-Cre':               ('vip',          ['cre']),
    'Pvalb-T2A-FlpO':             ('pvalb',        ['flp']),
    'Sst-IRES-FlpO':              ('sst',          ['flp']),
    'Vip-IRES-FlpO':              ('vip',          ['flp']),
    'Slc17a8-IRES2-Cre':          ('slc17a8',      ['cre']),
}

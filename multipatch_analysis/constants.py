from .genotypes import FLUOROPHORES

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

# note: the lower-case af488 and cascade_blue are for backward compatibility; these
# may be removed at some point
ALL_LABELS = ['biocytin', 'af488', 'cascade_blue'] + list(FLUOROPHORES.keys())

LAYERS = ['1', '2', '2/3', '3', '4', '5', '5a', '5b', '6']

ACSF_RECIPES = ["2mM Ca & Mg", "1.3mM Ca & 1mM Mg"]
INTERNAL_RECIPES = ["Standard K-Gluc", "PatchSeq", "K-Gluc -EGTA"]
INTERNAL_DYES = ['Cascade Blue', 'AF488', 'no dye']

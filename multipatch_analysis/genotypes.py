from __future__ import print_function
import sys, itertools

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
# The Genotype class below implements the second approach.


# Background research on understanding genotype symbols:
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
# dgCre : functions like Cre with reduced activity; increased by administration of TMP
# TRE : tetracycline response element - series of TetO sequences that
#       require tTA or rtTA to increase downstream expression
#         - dox binds tTa, which prevents further binding to TetO
#         - dox binds rtTa, which allows further binding to TetO
#       this system is also used for amplification.
# tTA : tet-off  (activates TRE only in absence of doxycycline)
# rtTA: tet-on   (activates TRE only in presence of doxycycline)
# tTA2 : same behavior as tTA

EXPRESSION_FACTORS = ['cre', 'flp', 'dre', 'tTA']

   
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
    'Pvalb-2A-FlpO':              ('pvalb',        ['flp']),
    'Cux2-CreERT2':               ('cux2',         ['cre']),
    'Chrna2-Cre_OE25':            ('chrna2',       ['cre']),
    'Penk-IRES2-Cre-neo':         ('penk',         ['cre']),
    'Slc17a6-IRES2-FlpO':         ('slc17a6',      ['flp']),
    'Ctgf-T2A-dgCre':             ('ctgf',         ['cre']),
}


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
    'Ai78(TITL-VSFPB)':           ('cre&tTA',    ['VSFP','Butterfly 1.2']),
    'Ai79(TITL-Jaws)':            ('cre&tTA',    ['Jaws','GFP','ER2']),
    'Ai82(TITL-GFP)':             ('cre&tTA',    ['EGFP']),
    'Ai85(TITL-iGluSnFR)':        ('cre&tTA',    ['iGluSnFR']),
    'Ai87(RCL-iGluSnFR)':         ('cre',        ['iGluSnFR']),
    'Ai92(TITL-YCX2.60)':         ('cre&tTA',    ['YCX2.60']),
    'Ai93(TITL-GCaMP6f)':         ('cre&tTA',    ['GCaMP6f']),
    'Ai94(TITL-GCaMP6s)':         ('cre&tTA',    ['GCaMP6s']),
    'Ai95(RCL-GCaMP6f)':          ('cre',        ['GCaMP6f']),
    'Ai96(RCL-GCaMP6s)':          ('cre',        ['GCaMP6s']),
    'Ai139(TIT2L-GFP-ICL-TPT)-D': ('cre',        ['EGFP','tdTomato']),
    'Ai139(TIT2L-GFP-ICL-TPT)':   ('cre',        ['EGFP','tdTomato']),
    'Ai140(TIT2L-GFP-ICL-tTA2)':  ('cre',        ['EGFP', 'tTA']),
    'Snap25-LSL-F2A-GFP':         ('cre',        ['EGFP']),
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


class Genotype(object):
    """
    Class used to provide phenotype information from a genotype string.

    Parameters
    ----------
    gtype : str
        A genotype string like
        'Tlx3-Cre_PL56/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt'

    Notes
    -----

    In the example given above, the genotype specifies two driver lines (Tlx3-Cre_PL65
    and Sst-IRES-FlpO) and two reporter lines (Ai65F and Ai140). This results in
    tdTomato expression in sst cells and GFP expression in tlx3 cells.

    Example usage::

        gtype_str = 'Tlx3-Cre_PL56/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt'
        gtype = Genotype(gtype_str)

        # Parsed genotype parts
        gtype.driver_lines                       # =>  ['Tlx3-Cre_PL65', 'Sst-IRES-FlpO']
        gtype.reporter_lines                     # =>  ['Ai65F', 'Ai140']

        # Reduced representation of genotype parts
        gtype.all_drivers                        # =>  ['tlx3', 'sst']
        gtype.all_reporters                      # =>  ['tdTomato', 'GFP']
        gtype.all_colors                         # =>  ['red', 'green']

        # What reporters are expressed if the cell is sst positive?
        gtype.expressed_reporters('tlx3')        # =>  ['tdTomato'] 

        # What colors are expressed if the cell is sst positive?
        gtype.expressed_colors('sst')            # =>  ['green']

        # Given colors seen in a cell, what combinations of drivers 
        # could have produced this?
        colors = {'red': True, 'green': False}
        gtype.predict_driver_expression(colors)  # =>  {('tlx3',): True, ('sst',): False, ('tlx3', 'sst'): False, (,): False}

    """
    def __init__(self, gtype):
        # describes which reporters are activated for each combination of drivers
        #    { (drivers,...): set([reporters_activated_by_drivers]), ... }
        self._driver_reporter_map = {}  
        self.gtype = gtype
        self._parse()

    def __repr__(self):
        return "<Genotype %s>" % self.gtype
    
    def expressed_reporters(self, drivers):
        """Return a set of reporters in this genotype (such as GFP, tdTomato, etc.)
        that would be expressed given a list of active drivers (such as pvalb, tlx3, etc.)

        Parameters
        ----------
        drivers : set
            The set of active drivers
        """
        return self._driver_reporter_map[tuple(sorted(drivers))]

    def expressed_colors(self, drivers):
        """Return a set of fluorescent emission colors generated by reporters in this genotype
        that would be expressed given a list of active drivers (such as pvalb, tlx3, etc.)

        Parameters
        ----------
        drivers : set
            The set of active drivers
        """
        reporters = self.expressed_reporters(drivers)
        return set([FLUOROPHORES[r] for r in reporters])

    def predict_driver_expression(self, colors):
        """Given information about fluorescent colors expressed in a cell,
        return predictions about whether each driver could have been active.

        Parameters
        ----------
        colors : dict
            Describes whether each color observed in a cell was expressed (True),
            not expressed (False), or ambiguous (None).

        Returns
        -------
        driver_expression : dict
            Dict indicating whether each driver in the genotype must be expressed
            (True), must not be expressed (False), or cannot be determined (None).

        Notes
        -----

        Example::

            colors = {
                'red':   True,   # cell is red
                'green': False,  # cell is not green
            }
        """
        driver_combos = self.test_driver_combinations(colors)
        true_combos = [drivers for drivers,prediction in driver_combos.items() if prediction is True]
        false_combos = [drivers for drivers,prediction in driver_combos.items() if prediction is False]
        driver_expression = {}
        for driver in self.all_drivers:
            driver_true_prediction_true = len([d for d in true_combos if driver in d])
            driver_false_prediction_true = len([d for d in true_combos if driver not in d])
            driver_true_prediction_false = len([d for d in false_combos if driver in d])
            driver_false_prediction_false = len([d for d in false_combos if driver not in d])

            if driver_true_prediction_true > 0 and driver_false_prediction_true == 0:
                # If every possible driver combination that is consistent with observed colors
                # contains this driver, then we asy it is definitely expressed.
                driver_expression[driver] = True
            elif driver_true_prediction_true == 0:
                # If no possible combinations that contain this driver are consistent with observed
                # colors, then we say the driver is definitely not expressed.
                driver_expression[driver] = False
            else:
                # Otherwise, we can't say one way or another whether this driver is expressed.
                driver_expression[driver] = None

        return driver_expression

    def test_driver_combinations(self, colors):
        """Given information about fluorescent colors expressed in a cell,
        return predictions about whether each combination of drivers could
        possibly have generated the specified color measurements.

        Parameters
        ----------
        colors : dict
            Describes whether each color observed in a cell was expressed (True),
            not expressed (False), or ambiguous (None). Example::

                colors = {
                    'red':   True,   # cell is red
                    'green': False,  # cell is not green
                    'blue':  None,   # cell may or may not be blue (no information, or ambiguous appearance)
                }
        
        Returns
        -------
        driver_expression : dict
            Dict indicating whether each driver combination in the genotype produces colors
            that are consistent (True) or inconsistent (False) with the observed colors.

        Notes
        -----

        The return values are only False where the predicted color expression is *inconsistent* with
        observed colors; in cases where there is not enough observed information, the return value
        will be True.
        """
        predictions = {}
        for drivers in self._driver_combinations():
            driver_combo_possible = True
            predicted_colors = self.expressed_colors(drivers=drivers)
            for color in self.all_colors:
                color_expressed = colors.get(color, None)
                if color_expressed is None:
                    continue
                color_predicted = color in predicted_colors
                if color_predicted != color_expressed:
                    driver_combo_possible = False
                    break
            predictions[drivers] = driver_combo_possible
        return predictions

    def _parse(self):
        """Attempt to predict phenotype information from a genotype string
        """
        ignore = ['wt', 'PhiC31-neo']
        
        parts = set()
        for part in self.gtype.split(';'):
            for subpart in part.split('/'):
                if subpart in ignore:
                    continue
                parts.add(subpart)
                
        self.driver_lines = [p for p in parts if p in DRIVER_LINES]
        self.reporter_lines = [p for p in parts if p in REPORTER_LINES]
        extra = parts - set(self.driver_lines + self.reporter_lines)
        if len(extra) > 0:
            raise Exception("Unknown genotype part(s): %s" % str(extra))

        self.all_drivers = set([DRIVER_LINES[d][0] for d in self.driver_lines])
        self.all_reporters = set()
        for r in self.reporter_lines:
            self.all_reporters |= set(REPORTER_LINES[r][1])
        self.all_reporters = self.all_reporters & set(FLUOROPHORES.keys())
        self.all_colors = set([FLUOROPHORES[r] for r in self.all_reporters])

        # generate a combined map describing input factors and the output products that would be generated
        #   e.g.:  [ (('tlx3',), ('cre')), (('cre',), ('tdTomato', 'tTA')) ]
        #             tlx3 generates cre,   cre generates tomato and tTA
        self._product_map = []
        for d in self.driver_lines:
            inputs, outputs = DRIVER_LINES[d]
            self._product_map.append((set(inputs.split('&')), set(outputs)))
        for r in self.reporter_lines:
            inputs, outputs = REPORTER_LINES[r]
            self._product_map.append((set(inputs.split('&')), set(outputs)))

        # Automatically determine driver-reporter mapping.
        self._driver_reporter_map = self._simulate_driver_combos()

        # generate the reverse mapping
        self._reporter_driver_map = {}
        for drivers, reporters in self._driver_reporter_map.items():
            for r in reporters:
                self._reporter_driver_map.setdefault(r, []).append(drivers)

    def _simulate_driver_combos(self):
        """Given all of the available driver lines in this genotype, simulate
        the reporter expression for each possible combination of drivers.

        Returns
        -------
        driver_reporter_map : dict
            Maps {(drivers,): [reporters]} to describe all of the reporters that would be expressed
            by each possible combination of drivers.
        """
        driver_reporter_map = {}

        # Generate a list of all possible driver line combinations
        #   e.g.:   [(), ('Tlx3-Cre_PL65',), ('Sst-IRES-FlpO',), ('Tlx3-Cre_PL65', 'Sst-IRES-FlpO')]
        driver_combos = self._driver_combinations()

        # Iterate over all combinations of the available driver lines and 
        # predict which reporters would be expressed in each case
        for drivers in driver_combos:
            expressed_factors = self._simulate_expression(drivers)
            expressed_reporters = set([f for f in expressed_factors if f in FLUOROPHORES])
            driver_reporter_map[tuple(sorted(drivers))] = expressed_reporters

        return driver_reporter_map

    def _driver_combinations(self):
        driver_line_combos = []
        for i in range(len(self.driver_lines) + 1):
            driver_line_combos.extend(list(itertools.combinations(self.driver_lines, i)))
        # convert driver line names to driver names  (e.g. 'Sst-IRES-FlpO' => 'sst')
        driver_combos = [tuple(sorted([DRIVER_LINES[d][0] for d in driver_lines])) for driver_lines in driver_line_combos]
        return driver_combos

    def _simulate_expression(self, factors):
        """Given a list of factors that could impact expression (drivers, recombinases, drugs, etc.),
        return a list of the reporters that would be expressed.

        Example:

            A pvalb-positive cell:  _simulate_expression(['pvalb'])
            An sst-positive cell plus tetracycline:  _simulate_expression(['sst', 'tet'])

        """
        factors = set(factors)
        new_expression = False
        for inputs, outputs in self._product_map:
            # inputs/outputs are tuples saying "if we have everything in inputs, then we generate everything in outputs"
            # e.g.:  inputs=('tlx3',)  outputs=('cre', 'tTA')
            #        inputs=('tTA',)   outputs=('tdTomato')
            if inputs.issubset(factors) and not outputs.issubset(factors):
                # Something new was expressed; add it to the list and run again
                factors |= outputs
                new_expression = True
        
        if new_expression:
            return self._simulate_expression(factors)
        else:
            return factors

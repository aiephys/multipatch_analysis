from multipatch_analysis import constants
from multipatch_analysis.genotypes import Genotype


def test_double_trans():
    # simple and easy: sst -> tomato
    gt = Genotype('Sst-IRES-Cre/wt;Ai14(RCL-tdT)/wt')
    assert set(gt.driver_lines) == set(['Sst-IRES-Cre'])
    assert set(gt.reporter_lines) == set(['Ai14(RCL-tdT)'])
    assert gt.all_drivers == set(['sst'])
    assert gt.all_reporters == set(['tdTomato'])
    assert gt.all_colors == set(['red'])

    assert gt.expressed_reporters([]) == set([])
    assert gt.expressed_reporters(['sst']) == set(['tdTomato'])
    
    assert gt.expressed_reporters([]) == set([])
    assert gt.expressed_reporters(['sst']) == set(['tdTomato'])
    
    assert gt.test_driver_combinations({'red': None}) == {(): True, ('sst',): True}
    assert gt.test_driver_combinations({'red': True}) == {(): False, ('sst',): True}
    assert gt.test_driver_combinations({'red': False}) == {(): True, ('sst',): False}

    assert gt.predict_driver_expression({'red': None}) == {'sst': None}
    assert gt.predict_driver_expression({'red': True}) == {'sst': True}
    assert gt.predict_driver_expression({'red': False}) == {'sst': False}


def test_quad_trans():
    # A little harder: tlx3 -> egfp, pvalb -> tomato
    gt = Genotype('Tlx3-Cre_PL56/wt;Pvalb-2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt')    
    assert set(gt.driver_lines) == set(['Tlx3-Cre_PL56', 'Pvalb-2A-FlpO'])
    assert set(gt.reporter_lines) == set(['Ai140(TIT2L-GFP-ICL-tTA2)', 'Ai65F'])
    assert gt.all_drivers == set(['tlx3', 'pvalb'])
    assert gt.all_reporters == set(['tdTomato', 'EGFP'])
    assert gt.all_colors == set(['green', 'red'])

    assert gt.expressed_reporters([]) == set([])
    assert gt.expressed_reporters(['pvalb']) == set(['tdTomato'])
    assert gt.expressed_reporters(['tlx3']) == set(['EGFP'])
    assert gt.expressed_reporters(['pvalb', 'tlx3']) == set(['tdTomato', 'EGFP'])

    assert gt.expressed_colors([]) == set()
    assert gt.expressed_colors(['pvalb']) == set(['red'])
    assert gt.expressed_colors(['tlx3']) == set(['green'])
    assert gt.expressed_colors(['pvalb', 'tlx3']) == set(['red', 'green'])

    assert gt.test_driver_combinations({}) == {(): True, ('pvalb',): True, ('tlx3',): True, ('pvalb', 'tlx3'): True}
    assert gt.test_driver_combinations({'red': True}) == {(): False, ('pvalb',): True, ('tlx3',): False, ('pvalb', 'tlx3'): True}
    assert gt.test_driver_combinations({'green': True}) == {(): False, ('pvalb',): False, ('tlx3',): True, ('pvalb', 'tlx3'): True}
    assert gt.test_driver_combinations({'red': False}) == {(): True, ('pvalb',): False, ('tlx3',): True, ('pvalb', 'tlx3'): False}
    assert gt.test_driver_combinations({'green': False}) == {(): True, ('pvalb',): True, ('tlx3',): False, ('pvalb', 'tlx3'): False}
    assert gt.test_driver_combinations({'green': False, 'red': True}) == {(): False, ('pvalb',): True, ('tlx3',): False, ('pvalb', 'tlx3'): False}
    assert gt.test_driver_combinations({'green': True, 'red': False}) == {(): False, ('pvalb',): False, ('tlx3',): True, ('pvalb', 'tlx3'): False}
    assert gt.test_driver_combinations({'green': False, 'red': False}) == {(): True, ('pvalb',): False, ('tlx3',): False, ('pvalb', 'tlx3'): False}
    assert gt.test_driver_combinations({'green': True, 'red': True}) == {(): False, ('pvalb',): False, ('tlx3',): False, ('pvalb', 'tlx3'): True}

    assert gt.predict_driver_expression({}) == {'tlx3': None, 'pvalb': None}
    assert gt.predict_driver_expression({'red': True}) == {'tlx3': None, 'pvalb': True} 
    assert gt.predict_driver_expression({'green': True}) == {'tlx3': True, 'pvalb': None} 
    assert gt.predict_driver_expression({'red': False}) == {'tlx3': None, 'pvalb': False} 
    assert gt.predict_driver_expression({'green': False}) == {'tlx3': False, 'pvalb': None} 
    assert gt.predict_driver_expression({'green': False, 'red': True}) == {'tlx3': False, 'pvalb': True} 
    assert gt.predict_driver_expression({'green': True, 'red': False}) == {'tlx3': True, 'pvalb': False} 
    assert gt.predict_driver_expression({'green': False, 'red': False}) == {'tlx3': False, 'pvalb': False} 
    assert gt.predict_driver_expression({'green': True, 'red': True}) == {'tlx3': True, 'pvalb': True} 


def test_overlapping_quad():
    # harder yet: pv -> egfp+tomato, rorb -> tomato
    gt = Genotype('Pvalb-IRES-Cre/wt;Rorb-T2A-tTA2/wt;Ai63(TIT-tdT)/Ai140(TIT2L-GFP-ICL-tTA2)')    
    assert set(gt.driver_lines) == set(['Rorb-T2A-tTA2', 'Pvalb-IRES-Cre'])
    assert set(gt.reporter_lines) == set(['Ai140(TIT2L-GFP-ICL-tTA2)', 'Ai63(TIT-tdT)'])
    assert gt.all_drivers == set(['rorb', 'pvalb'])
    assert gt.all_reporters == set(['tdTomato', 'EGFP'])
    assert gt.all_colors == set(['green', 'red'])

    assert gt.expressed_reporters([]) == set([])
    assert gt.expressed_reporters(['rorb']) == set(['tdTomato'])
    assert set(gt.expressed_reporters(['pvalb'])) == set(['tdTomato', 'EGFP'])
    assert gt.expressed_reporters(['rorb', 'pvalb']) == set(['tdTomato', 'EGFP'])
    assert gt.expressed_reporters(['pvalb', 'rorb']) == set(['tdTomato', 'EGFP'])

    assert gt.test_driver_combinations({}) == {(): True, ('rorb',): True, ('pvalb',): True, ('pvalb', 'rorb'): True}
    assert gt.test_driver_combinations({'red': True}) == {(): False, ('rorb',): True, ('pvalb',): True, ('pvalb', 'rorb'): True}
    assert gt.test_driver_combinations({'green': True}) == {(): False, ('rorb',): False, ('pvalb',): True, ('pvalb', 'rorb'): True}
    assert gt.test_driver_combinations({'red': False}) == {(): True, ('rorb',): False, ('pvalb',): False, ('pvalb', 'rorb'): False}
    assert gt.test_driver_combinations({'green': False}) == {(): True, ('rorb',): True, ('pvalb',): False, ('pvalb', 'rorb'): False}
    assert gt.test_driver_combinations({'green': False, 'red': True}) == {(): False, ('rorb',): True, ('pvalb',): False, ('pvalb', 'rorb'): False}
    assert gt.test_driver_combinations({'green': True, 'red': False}) == {(): False, ('rorb',): False, ('pvalb',): False, ('pvalb', 'rorb'): False}
    assert gt.test_driver_combinations({'green': False, 'red': False}) == {(): True, ('rorb',): False, ('pvalb',): False, ('pvalb', 'rorb'): False}
    assert gt.test_driver_combinations({'green': True, 'red': True}) == {(): False, ('rorb',): False, ('pvalb',): True, ('pvalb', 'rorb'): True}

    assert gt.predict_driver_expression({}) == {'pvalb': None, 'rorb': None}
    assert gt.predict_driver_expression({'red': True}) == {'pvalb': None, 'rorb': None} 
    assert gt.predict_driver_expression({'green': True}) == {'pvalb': True, 'rorb': None} 
    assert gt.predict_driver_expression({'red': False}) == {'pvalb': False, 'rorb': False} 
    assert gt.predict_driver_expression({'green': False}) == {'pvalb': False, 'rorb': None} 
    assert gt.predict_driver_expression({'green': False, 'red': True}) == {'pvalb': False, 'rorb': True} 
    assert gt.predict_driver_expression({'green': True, 'red': False}) == {'pvalb': False, 'rorb': False} 
    assert gt.predict_driver_expression({'green': False, 'red': False}) == {'pvalb': False, 'rorb': False} 
    assert gt.predict_driver_expression({'green': True, 'red': True}) == {'pvalb': True, 'rorb': None} 

    assert gt.expressed_colors([]) == set()
    assert gt.expressed_colors(['rorb']) == set(['red'])
    assert gt.expressed_colors(['pvalb']) == set(['red', 'green'])
    assert gt.expressed_colors(['pvalb', 'rorb']) == set(['red', 'green'])


def test_known_genotypes():
    # make sure simulation matches known results
    for gtyp_str, mapping in constants.GENOTYPES.items():
        gt = Genotype(gtyp_str)
        driver_reporter_map = gt._simulate_driver_combos()
        for k,v in mapping.items():
            assert set(driver_reporter_map[k]) == set(v)



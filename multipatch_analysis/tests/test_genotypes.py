from multipatch_analysis import constants
from multipatch_analysis.genotypes import Genotype


def test_genotypes():
    gt = Genotype('Pvalb-IRES-Cre/wt;Rorb-T2A-tTA2/wt;Ai63(TIT-tdT)/Ai140(TIT2L-GFP-ICL-tTA2)')
    assert set(gt.driver_lines) == set(['Rorb-T2A-tTA2', 'Pvalb-IRES-Cre'])
    assert set(gt.reporter_lines) == set(['Ai140(TIT2L-GFP-ICL-tTA2)', 'Ai63(TIT-tdT)'])

    assert gt.all_drivers == set(['rorb', 'pvalb'])
    assert gt.all_reporters == set(['tdTomato', 'EGFP'])
    assert gt.all_colors == set(['green', 'red'])

    assert gt.expressed_reporters([]) == []
    assert gt.expressed_reporters(['rorb']) == ['tdTomato']
    assert set(gt.expressed_reporters(['pvalb'])) == set(['tdTomato', 'EGFP'])
    assert gt.expressed_reporters(['rorb', 'pvalb']) == ['tdTomato', 'EGFP']
    assert gt.expressed_reporters(['pvalb', 'rorb']) == ['tdTomato', 'EGFP']

    assert gt.predict_driver_expression({}) == {(): None, ('rorb',): None, ('pvalb',): None, ('pvalb', 'rorb'): None}
    assert gt.predict_driver_expression({'red': True}) == {(): False, ('rorb',): True, ('pvalb',): True, ('pvalb', 'rorb'): True}
    assert gt.predict_driver_expression({'green': True}) == {(): False, ('rorb',): False, ('pvalb',): True, ('pvalb', 'rorb'): True}
    assert gt.predict_driver_expression({'red': False}) == {(): None, ('rorb',): False, ('pvalb',): False, ('pvalb', 'rorb'): False}
    assert gt.predict_driver_expression({'green': False}) == {(): None, ('rorb',): None, ('pvalb',): False, ('pvalb', 'rorb'): False}
    assert gt.predict_driver_expression({'green': False, 'red': True}) == {(): False, ('rorb',): True, ('pvalb',): False, ('pvalb', 'rorb'): False}
    assert gt.predict_driver_expression({'green': True, 'red': False}) == {(): False, ('rorb',): False, ('pvalb',): False, ('pvalb', 'rorb'): False}
    assert gt.predict_driver_expression({'green': False, 'red': False}) == {(): None, ('rorb',): False, ('pvalb',): False, ('pvalb', 'rorb'): False}
    assert gt.predict_driver_expression({'green': True, 'red': True}) == {(): False, ('rorb',): False, ('pvalb',): True, ('pvalb', 'rorb'): True}

    assert gt.expressed_colors([]) == set()
    assert gt.expressed_colors(['rorb']) == set(['red'])
    assert gt.expressed_colors(['pvalb']) == set(['red', 'green'])

    for gtyp_str, mapping in constants.GENOTYPES.items():
        gt = Genotype(gtyp_str)
        driver_reporter_map = gt._simulate_driver_combos()
        for k,v in mapping.items():
            assert set(driver_reporter_map[k]) == set(v)



from aisynphys.genotypes import Genotype


known_genotypes = dict([
    ('Chat-IRES-Cre-neo/wt;Snap25-LSL-F2A-GFP/wt', {(): [], ('chat',): ['EGFP']}),
    ('Chrna2-Cre_OE25/wt;Ai14(RCL-tdT)/wt', {(): [], ('chrna2',): ['tdTomato']}),
    ('Ctgf-T2A-dgCre/wt;Ai14(RCL-tdT)/wt', {('ctgf',): ['tdTomato'], (): []}),
    ('Cux2-CreERT2/wt;Ai14(RCL-tdT)/wt', {('cux2',): ['tdTomato'], (): []}),
    ('Nr5a1-Cre/wt;Ai14(RCL-tdT)/wt', {(): [], ('nr5a1',): ['tdTomato']}),
    ('Ntsr1-Cre_GN220/wt;Ai14(RCL-tdT)/wt', {('ntsr1',): ['tdTomato'], (): []}),
    ('Ntsr1-Cre_GN220/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt;Ai65F/wt', {('ntsr1',): ['EGFP'], (): []}),
    ('Ntsr1-Cre_GN220/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('ntsr1',): ['EGFP'], (): []}),
    ('Ntsr1-Cre_GN220/wt;Pvalb-T2A-FlpO/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('ntsr1',): ['EGFP'], (): [], ('ntsr1', 'pvalb'): ['EGFP'], ('pvalb',): []}),
    ('Ntsr1-Cre_GN220/wt;Pvalb-T2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('ntsr1',): ['EGFP'], (): [], ('ntsr1', 'pvalb'): ['tdTomato', 'EGFP'], ('pvalb',): ['tdTomato']}),
    ('Ntsr1-Cre_GN220/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('ntsr1',): ['EGFP'], (): [], ('ntsr1', 'pvalb'): ['tdTomato', 'EGFP'], ('pvalb',): ['tdTomato']}),
    ('Ntsr1-Cre_GN220/wt;Sst-IRES-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('ntsr1',): ['EGFP'], (): [], ('sst',): ['tdTomato'], ('ntsr1', 'sst'): ['tdTomato', 'EGFP']}),
    ('Ntsr1-Cre_GN220/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('ntsr1',): ['EGFP'], (): [], ('sst',): ['tdTomato'], ('ntsr1', 'sst'): ['tdTomato', 'EGFP']}),
    ('Ntsr1-Cre_GN220/wt;Vip-IRES-FlpO/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt;Ai65F/wt', {('ntsr1',): ['EGFP'], ('vip',): ['tdTomato'], (): [], ('ntsr1', 'vip'): ['tdTomato', 'EGFP']}),
    ('Ntsr1-Cre_GN220/wt;Vip-IRES-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('ntsr1',): ['EGFP'], ('vip',): ['tdTomato'], (): [], ('ntsr1', 'vip'): ['tdTomato', 'EGFP']}),
    ('Ntsr1-Cre_GN220/wt;Vip-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('ntsr1',): ['EGFP'], ('vip',): ['tdTomato'], (): [], ('ntsr1', 'vip'): ['tdTomato', 'EGFP']}),
    ('Penk-IRES2-Cre-neo/wt;Slc17a6-IRES2-FlpO/wt;Ai65(RCFL-tdT)/wt', {(): [], ('penk',): [], ('slc17a6',): [], ('penk', 'slc17a6'): ['tdTomato']}),
    ('Pvalb-IRES-Cre/wt;Ai14(RCL-tdT)/wt', {('pvalb',): ['tdTomato'], (): []}),
    ('Pvalb-IRES-Cre/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['EGFP'], (): []}),
    ('Pvalb-IRES-Cre/wt;Ai63(TIT-tdT)/Ai140(TIT2L-GFP-ICL-tTA2)', {('pvalb',): ['EGFP', 'tdTomato'], (): []}),
    ('Pvalb-IRES-Cre/wt;Rorb-T2A-tTA2/wt;Ai63(TIT-tdT)/Ai140(TIT2L-GFP-ICL-tTA2)', {('pvalb',): ['EGFP', 'tdTomato'], ('rorb',): ['tdTomato'], (): [], ('pvalb', 'rorb'): ['tdTomato', 'EGFP']}),
    ('Pvalb-IRES-Cre/wt;Rorb-T2A-tTA2/wt;Ai63(TIT-tdT)/wt', {('pvalb',): [], ('rorb',): ['tdTomato'], (): [], ('pvalb', 'rorb'): ['tdTomato']}),
    ('Pvalb-T2A-FlpO/wt;Ai65F/Ai65F', {('pvalb',): ['tdTomato'], (): []}),
    ('Pvalb-T2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], (): []}),
    ('Pvalb-T2A-FlpO/wt;Ai65F/wt', {('pvalb',): ['tdTomato'], (): []}),
    ('Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], (): []}),
    ('Rbp4-Cre_KL100/wt;Ai14(RCL-tdT)/wt', {(): [], ('rbp4',): ['tdTomato']}),
    ('Rorb-T2A-tTA2/wt;Ai63(TIT-tdT)/Ai140(TIT2L-GFP-ICL-tTA2)', {('rorb',): ['tdTomato'], (): []}),
    ('Sim1-Cre_KJ18/wt;Ai14(RCL-tdT)/wt', {('sim1',): ['tdTomato'], (): []}),
    ('Sim1-Cre_KJ18/wt;Ai65F/wt;Ai139(TIT2L-GFP-ICL-TPT)/wt', {('sim1',): ['tdTomato', 'EGFP'], (): []}),
    ('Sim1-Cre_KJ18/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('sim1',): ['EGFP'], (): [], ('pvalb', 'sim1'): ['tdTomato', 'EGFP']}),
    ('Sim1-Cre_KJ18/wt;Pvalb-T2A-FlpO/wt;Snap25-LSL-F2A-GFP/wt;Ai65F/wt', {('pvalb',): ['tdTomato'], ('sim1',): ['EGFP'], (): [], ('pvalb', 'sim1'): ['tdTomato', 'EGFP']}),
    ('Sim1-Cre_KJ18/wt;Sst-IRES-FlpO/wt;Ai139(TIT2L-GFP-ICL-TPT)/wt', {('sim1',): ['tdTomato', 'EGFP'], (): [], ('sim1', 'sst'): ['tdTomato', 'EGFP'], ('sst',): []}),
    ('Sim1-Cre_KJ18/wt;Sst-IRES-FlpO/wt;Ai65F/Ai65F;Ai139(TIT2L-GFP-ICL-TPT)/wt', {('sim1',): ['tdTomato', 'EGFP'], (): [], ('sim1', 'sst'): ['tdTomato', 'EGFP'], ('sst',): ['tdTomato']}),
    ('Sim1-Cre_KJ18/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai139(TIT2L-GFP-ICL-TPT)/wt', {('sim1',): ['tdTomato', 'EGFP'], (): [], ('sim1', 'sst'): ['tdTomato', 'EGFP'], ('sst',): ['tdTomato']}),
    ('Sim1-Cre_KJ18/wt;Sst-IRES-FlpO/wt;PhiC31-neo/Ai65F', {('sim1',): [], (): [], ('sim1', 'sst'): ['tdTomato'], ('sst',): ['tdTomato']}),
    ('Slc17a8-IRES2-Cre/wt;Ai14(RCL-tdT)/wt', {(): [], ('slc17a8',): ['tdTomato']}),
    ('Sst-IRES-Cre/wt;Ai14(RCL-tdT)/wt', {(): [], ('sst',): ['tdTomato']}),
    ('Sst-IRES-Cre/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {(): [], ('sst',): ['EGFP']}),
    ('Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/Pvalb-T2A-FlpO;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('pvalb', 'sst'): ['tdTomato', 'EGFP'], (): [], ('sst',): ['EGFP']}),
    ('Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/Pvalb-T2A-FlpO;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('pvalb', 'sst'): ['tdTomato', 'EGFP'], (): [], ('sst',): ['EGFP']}),
    ('Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): [], ('pvalb', 'sst'): ['EGFP'], (): [], ('sst',): ['EGFP']}),
    ('Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/Ai65F', {('pvalb',): ['tdTomato'], ('pvalb', 'sst'): ['tdTomato'], (): [], ('sst',): []}),
    ('Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('pvalb', 'sst'): ['tdTomato', 'EGFP'], (): [], ('sst',): ['EGFP']}),
    ('Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt', {('pvalb',): ['tdTomato'], ('pvalb', 'sst'): ['tdTomato'], (): [], ('sst',): []}),
    ('Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/Ai140(TIT2L-GFP-ICL-tTA2)', {('pvalb',): ['tdTomato'], ('pvalb', 'sst'): ['tdTomato', 'EGFP'], (): [], ('sst',): ['EGFP']}),
    ('Sst-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('pvalb', 'sst'): ['tdTomato', 'EGFP'], (): [], ('sst',): ['EGFP']}),
    ('Sst-IRES-FlpO/wt;Ai65F/Ai65F', {(): [], ('sst',): ['tdTomato']}),
    ('Sst-IRES-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {(): [], ('sst',): ['tdTomato']}),
    ('Sst-IRES-FlpO/wt;Ai65F/wt', {(): [], ('sst',): ['tdTomato']}),
    ('Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {(): [], ('sst',): ['tdTomato']}),
    ('Sst-IRES-FlpO/wt;PhiC31-neo/Ai65F', {(): [], ('sst',): ['tdTomato']}),
    ('Tlx3-Cre_PL56/wt;Ai14(RCL-tdT)/wt', {('tlx3',): ['tdTomato'], (): []}),
    ('Tlx3-Cre_PL56/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('tlx3',): ['EGFP'], (): []}),
    ('Tlx3-Cre_PL56/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('tlx3',): ['EGFP'], (): []}),
    ('Tlx3-Cre_PL56/wt;Pvalb-2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('tlx3',): ['EGFP'], (): [], ('pvalb', 'tlx3'): ['tdTomato', 'EGFP']}),
    ('Tlx3-Cre_PL56/wt;Pvalb-2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('tlx3',): ['EGFP'], (): [], ('pvalb', 'tlx3'): ['tdTomato', 'EGFP']}),
    ('Tlx3-Cre_PL56/wt;Pvalb-T2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('tlx3',): ['EGFP'], (): [], ('pvalb', 'tlx3'): ['tdTomato', 'EGFP']}),
    ('Tlx3-Cre_PL56/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt', {('pvalb',): ['tdTomato'], ('tlx3',): [], (): [], ('pvalb', 'tlx3'): ['tdTomato']}),
    ('Tlx3-Cre_PL56/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('tlx3',): ['EGFP'], (): [], ('pvalb', 'tlx3'): ['tdTomato', 'EGFP']}),
    ('Tlx3-Cre_PL56/wt;Pvalb-T2A-FlpO/wt;PhiC31-neo/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('tlx3',): ['EGFP'], (): [], ('pvalb', 'tlx3'): ['tdTomato', 'EGFP']}),
    ('Tlx3-Cre_PL56/wt;Sst-IRES-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('tlx3',): ['EGFP'], (): [], ('sst',): ['tdTomato'], ('sst', 'tlx3'): ['tdTomato', 'EGFP']}),
    ('Tlx3-Cre_PL56/wt;Sst-IRES-FlpO/wt;Ai65F/PhiC31-neo;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('tlx3',): ['EGFP'], (): [], ('sst',): ['tdTomato'], ('sst', 'tlx3'): ['tdTomato', 'EGFP']}),
    ('Tlx3-Cre_PL56/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('tlx3',): ['EGFP'], (): [], ('sst',): ['tdTomato'], ('sst', 'tlx3'): ['tdTomato', 'EGFP']}),
    ('Tlx3-Cre_PL56/wt;Sst-IRES-FlpO/wt;PhiC31-neo/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('tlx3',): ['EGFP'], (): [], ('sst',): ['tdTomato'], ('sst', 'tlx3'): ['tdTomato', 'EGFP']}),
    ('Tlx3-Cre_PL56/wt;Vip-IRES-FlpO/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('tlx3',): ['EGFP'], ('vip',): [], (): [], ('tlx3', 'vip'): ['EGFP']}),
    ('Tlx3-Cre_PL56/wt;Vip-IRES-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('tlx3',): ['EGFP'], ('vip',): ['tdTomato'], (): [], ('tlx3', 'vip'): ['tdTomato', 'EGFP']}),
    ('Tlx3-Cre_PL56/wt;Vip-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('tlx3',): ['EGFP'], ('vip',): ['tdTomato'], (): [], ('tlx3', 'vip'): ['tdTomato', 'EGFP']}),
    ('Tlx3-Cre_PL56/wt;Vip-IRES-FlpO/wt;PhiC31-neo/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('tlx3',): ['EGFP'], ('vip',): ['tdTomato'], (): [], ('tlx3', 'vip'): ['tdTomato', 'EGFP']}),
    ('Vip-IRES-Cre/wt;Ai14(RCL-tdT)/wt', {('vip',): ['tdTomato'], (): []}),
    ('Vip-IRES-Cre/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('vip',): ['EGFP'], (): []}),
    ('Vip-IRES-Cre/wt;Pvalb-T2A-FlpO/Pvalb-T2A-FlpO;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('vip',): ['EGFP'], (): [], ('pvalb', 'vip'): ['tdTomato', 'EGFP']}),
    ('Vip-IRES-Cre/wt;Pvalb-T2A-FlpO/Pvalb-T2A-FlpO;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('vip',): ['EGFP'], (): [], ('pvalb', 'vip'): ['tdTomato', 'EGFP']}),
    ('Vip-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('vip',): ['EGFP'], (): [], ('pvalb', 'vip'): ['tdTomato', 'EGFP']}),
    ('Vip-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/Ai140(TIT2L-GFP-ICL-tTA2)', {('pvalb',): ['tdTomato'], ('vip',): ['EGFP'], (): [], ('pvalb', 'vip'): ['tdTomato', 'EGFP']}),
    ('Vip-IRES-Cre/wt;Pvalb-T2A-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('pvalb',): ['tdTomato'], ('vip',): ['EGFP'], (): [], ('pvalb', 'vip'): ['tdTomato', 'EGFP']}),
    ('Vip-IRES-Cre/wt;Sst-IRES-FlpO/wt;Ai65F/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('sst', 'vip'): ['tdTomato', 'EGFP'], ('vip',): ['EGFP'], (): [], ('sst',): ['tdTomato']}),
    ('Vip-IRES-Cre/wt;Sst-IRES-FlpO/wt;Ai65F/wt;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('sst', 'vip'): ['tdTomato', 'EGFP'], ('vip',): ['EGFP'], (): [], ('sst',): ['tdTomato']}),
    ('Vip-IRES-Cre/wt;Sst-IRES-FlpO/wt;PhiC31-neo/Ai65F;Ai140(TIT2L-GFP-ICL-tTA2)/wt', {('sst', 'vip'): ['tdTomato', 'EGFP'], ('vip',): ['EGFP'], (): [], ('sst',): ['tdTomato']}),

])


def test_known_genotypes():
    # make sure simulation matches known results
    for gtyp_str, mapping in known_genotypes.items():
        gt = Genotype(gtyp_str)
        for drivers, reporters in mapping.items():
            assert gt.expressed_reporters(drivers) == set(reporters)


def test_single_trans():
    gt = Genotype('Ai2/wt')

    assert gt.driver_lines == []
    assert gt.reporter_lines == ['Ai2']
    assert gt.all_drivers == set()
    assert gt.all_reporters == set(['EYFP'])
    assert gt.all_colors == set(['yellow'])

    assert gt.expressed_reporters([]) == set(['EYFP'])
    assert gt.expressed_colors([]) == set(['yellow'])
    assert gt.predict_driver_expression({'yellow': True}) == {}



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
    
    assert gt.expressed_colors([]) == set([])
    assert gt.expressed_colors(['sst']) == set(['red'])
    
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
    assert gt.expressed_reporters(['pvalb']) == set(['tdTomato', 'EGFP'])
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


def test_intersectional():
    # penk+slc17a6 => tomato
    gt = Genotype('Penk-IRES2-Cre-neo/wt;Slc17a6-IRES2-FlpO/wt;Ai65(RCFL-tdT)/wt')
    assert set(gt.driver_lines) == set(['Penk-IRES2-Cre-neo', 'Slc17a6-IRES2-FlpO'])
    assert set(gt.reporter_lines) == set(['Ai65(RCFL-tdT)'])
    assert gt.all_drivers == set(['penk', 'slc17a6'])
    assert gt.all_reporters == set(['tdTomato'])
    assert gt.all_colors == set(['red'])

    assert gt.expressed_reporters([]) == set([])
    assert gt.expressed_reporters(['slc17a6']) == set([])
    assert gt.expressed_reporters(['penk']) == set([])
    assert gt.expressed_reporters(['slc17a6', 'penk']) == set(['tdTomato'])
    assert gt.expressed_reporters(['penk', 'slc17a6']) == set(['tdTomato'])

    assert gt.expressed_colors([]) == set([])
    assert gt.expressed_colors(['slc17a6']) == set([])
    assert gt.expressed_colors(['penk']) == set([])
    assert gt.expressed_colors(['slc17a6', 'penk']) == set(['red'])
    assert gt.expressed_colors(['penk', 'slc17a6']) == set(['red'])

    assert gt.test_driver_combinations({}) == {(): True, ('slc17a6',): True, ('penk',): True, ('penk', 'slc17a6'): True}
    assert gt.test_driver_combinations({'red': None}) == {(): True, ('slc17a6',): True, ('penk',): True, ('penk', 'slc17a6'): True}
    assert gt.test_driver_combinations({'red': True}) == {(): False, ('slc17a6',): False, ('penk',): False, ('penk', 'slc17a6'): True}
    assert gt.test_driver_combinations({'red': False}) == {(): True, ('slc17a6',): True, ('penk',): True, ('penk', 'slc17a6'): False}

    assert gt.predict_driver_expression({}) == {'penk': None, 'slc17a6': None}
    assert gt.predict_driver_expression({'red': True}) == {'penk': True, 'slc17a6': True} 
    assert gt.predict_driver_expression({'red': False}) == {'penk': None, 'slc17a6': None} 


def test_dox():
    gt = Genotype('Ai63(TIT-tdT)/wt;Rorb-T2A-tTA2/wt')
    assert set(gt.driver_lines) == set(['Rorb-T2A-tTA2'])
    assert set(gt.reporter_lines) == set(['Ai63(TIT-tdT)'])
    assert gt.all_drivers == set(['rorb'])
    assert gt.all_reporters == set(['tdTomato'])
    assert gt.all_colors == set(['red'])

    assert gt.expressed_reporters([]) == set([])
    assert gt.expressed_reporters(['rorb']) == set(['tdTomato'])
    assert gt.expressed_reporters(['rorb', 'dox']) == set([])

    assert gt.expressed_colors([]) == set([])
    assert gt.expressed_colors(['rorb']) == set(['red'])
    assert gt.expressed_colors(['rorb', 'dox']) == set([])

    assert gt.test_driver_combinations({}) == {(): True, ('rorb',): True}
    assert gt.test_driver_combinations({'red': True}) == {(): False, ('rorb',): True}
    assert gt.test_driver_combinations({'red': False}) == {(): True, ('rorb',): False}
    assert gt.test_driver_combinations({'red': True}, starting_factors=['dox']) == {(): False, ('rorb',): False}
    assert gt.test_driver_combinations({'red': False}, starting_factors=['dox']) == {(): True, ('rorb',): True}

    assert gt.predict_driver_expression({}) == {'rorb': None}
    assert gt.predict_driver_expression({'red': True}) == {'rorb': True}
    assert gt.predict_driver_expression({'red': False}) == {'rorb': False}
    # TODO: Logic engine can't handle this case yet!
    # assert gt.predict_driver_expression({'red': True}, starting_factors=['dox']) == {'rorb': None}
    assert gt.predict_driver_expression({'red': False}, starting_factors=['dox']) == {'rorb': None}


from pytest import raises
from aisynphys.cell_class import CellClass
from aisynphys.database import default_db as db

def test_cell_class_init():
    with raises(AssertionError):
        CellClass(asdf=1)

    with raises(AssertionError):
        CellClass(db.Pair.n_ex_test_spikes < 10)

    with raises(AssertionError):
        CellClass([])

    with raises(AssertionError):
        CellClass(name=1)


def test_cell_class_is_excitatory():

    cls = CellClass(cre_type='sst')
    assert cls.is_excitatory is False
    assert cls.output_synapse_type == 'in'

    cls = CellClass(cre_type='sim1')
    assert cls.is_excitatory is True
    assert cls.output_synapse_type == 'ex'

    cls = CellClass(cre_type=['sst', 'pvalb'])
    assert cls.is_excitatory is False
    assert cls.output_synapse_type == 'in'

    cls = CellClass(cre_type=['sst', 'tlx3'])
    assert cls.is_excitatory is None
    assert cls.output_synapse_type is None

    cls = CellClass(cre_type=['sst'], dendrite_type='spiny')
    assert cls.is_excitatory is None
    assert cls.output_synapse_type is None

    cls = CellClass(cre_type='rorb', dendrite_type='spiny')
    assert cls.is_excitatory is True
    assert cls.output_synapse_type == 'ex'

    cls = CellClass(cre_type='rorb', dendrite_type='sparsely spiny')
    assert cls.is_excitatory is None
    assert cls.output_synapse_type is None

    cls = CellClass(cre_type=['sst', 'vip'], dendrite_type='sparsely spiny')
    assert cls.is_excitatory is False
    assert cls.output_synapse_type == 'in'

    cls = CellClass(cell_class='in')
    assert cls.is_excitatory is False
    assert cls.output_synapse_type == 'in'

    cls = CellClass(cell_class_nonsynaptic='in')
    assert cls.is_excitatory is False
    assert cls.output_synapse_type == 'in'

    cls = CellClass(cell_class='ex')
    assert cls.is_excitatory is True
    assert cls.output_synapse_type == 'ex'

    cls = CellClass(cell_class_nonsynaptic='ex')
    assert cls.is_excitatory is True
    assert cls.output_synapse_type == 'ex'

    cls = CellClass(cell_class='in', cell_class_nonsynaptic='ex')
    assert cls.is_excitatory is None
    assert cls.output_synapse_type is None

    cls = CellClass(target_layer='2/3')
    assert cls.is_excitatory is None
    assert cls.output_synapse_type is None

    cls = CellClass(target_layer='2/3', cell_class='ex')
    assert cls.is_excitatory is True
    assert cls.output_synapse_type == 'ex'


def test_cell_class_contains():

    cell = db.Cell()
    morpho = db.Morphology()
    loc = db.CorticalCellLocation()
    cell.morphology = morpho
    cell.cortical_location = loc

    cls = CellClass(cre_type=('sst', 'pvalb'))
    
    cell.cre_type = 'sst'
    assert cell in cls

    cell.cre_type = 'sim1'
    assert cell not in cls

    cls = CellClass(cre_type='sim1', dendrite_type='spiny')

    morpho.dendrite_type = 'spiny'
    assert cell in cls
    
    morpho.dendrite_type = 'aspiny'
    assert cell not in cls
    
    cls = CellClass(cre_type='sim1', dendrite_type='spiny', cortical_layer='2/3')

    morpho.dendrite_type = 'spiny'
    loc.cortical_layer = '2/3'
    assert cell in cls

    loc.cortical_layer = '2'
    assert cell not in cls

    cls = CellClass(cre_type='sim1', dendrite_type='spiny', cortical_layer='2/3', cell_class='ex')

    cell.cell_class = 'ex'
    loc.cortical_layer = '2/3'
    assert cell in cls

    cell.cell_class = 'mixed'
    assert cell not in cls

    cls = CellClass(input_resistance=0)
    # no intrinsic added to this cell
    assert cell not in cls


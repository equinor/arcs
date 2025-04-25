import numpy as np

from arcs.model import _find_enclosing, run_reaction_calc, _get_gibbs_constant


def test_find_enclosing():
    values = [100, 200, 300, 400, 500]

    # Test with a value in between two existing values
    lower, upper = _find_enclosing(250, values)
    assert lower == 200
    assert upper == 300


def test_run_reaction_calc():
    temperature = 311
    pressure = 2

    # Assuming this function returns a list with the Gibbs free energy values
    gibbs_values, equilibrium, reaction_ids = run_reaction_calc(pressure, temperature)

    # Check if the returned values are of expected types
    assert isinstance(gibbs_values, list)
    assert isinstance(equilibrium, list)
    assert isinstance(reaction_ids, list)

    # Check if gibbs_values is not empty
    assert len(gibbs_values) > 0
    assert len(equilibrium) > 0
    assert len(reaction_ids) > 0

    # Check if the lengths of gibbs_values and reaction_ids match
    assert len(gibbs_values) == len(reaction_ids)


def test_get_gibbs_constant():
    # Mock data for reactions and pressure/temperature combinations
    reactions = [
        {1: {"e": None, "k": 1.0, "g": 1.0}, 2: {"e": None, "k": 1.5, "g": 10.0}},
        {1: {"e": None, "k": 1.0, "g": 2.0}, 2: {"e": None, "k": 1.5, "g": 20.0}},
        {1: {"e": None, "k": 1.0, "g": 3.0}, 2: {"e": None, "k": 1.5, "g": 30.0}},
        {1: {"e": None, "k": 2.0, "g": 4.0}, 2: {"e": None, "k": 2.5, "g": 40.0}},
    ]
    pt_combinations = [(350, 30), (350, 35), (400, 30), (400, 35)]
    pressure = 31
    temperature = 375

    # Call the function
    calculated_gibbs_constants, equilibrium, reaction_ids = _get_gibbs_constant(
        reactions, pt_combinations, pressure, temperature
    )

    # Check if the returned values are of expected types
    assert isinstance(calculated_gibbs_constants, list)
    assert isinstance(equilibrium, list)
    assert isinstance(reaction_ids, list)

    # Check if calculated_gibbs_constants matches expected length
    assert len(calculated_gibbs_constants) == len(reaction_ids)

    # Check if equilibrium has the expected number of entries
    assert len(equilibrium) > 0

    # Assert interpolated values are correct
    assert calculated_gibbs_constants == [np.array(2.2), np.array(22.)]

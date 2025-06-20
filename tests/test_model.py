import pytest
import numpy as np
from arcs.model import (
    TEMPERATURES,
    ReactionType,
    _find_enclosing,
    _interp_gibbs,
    interpolate_gibbs_values,
)


@pytest.mark.parametrize("value, expect_lower, expect_upper",[
    pytest.param(100, 100, 100, id="100 (lower bound exact match)"),
    pytest.param(500, 500, 500, id="500 (upper bound exact match)"),
    pytest.param(300, 300, 300, id="300 (exact match)"),
    pytest.param(250, 200, 300, id="250 (between 200 and 300)"),
    pytest.param(0, None, None, id="0 (out of bounds)"),
    pytest.param(600, None, None, id="600 (out of bounds)"),
])
def test_find_enclosing(value: int, expect_lower: int | None, expect_upper: int | None):
    values = np.array([100, 200, 300, 400, 500])

    if expect_lower is not None and expect_upper is not None:
        lower, upper = _find_enclosing(values, value)
        assert lower == expect_lower
        assert upper == expect_upper
    else:
        with pytest.raises(IndexError, match="Search value must be between 100 and 500"):
            _find_enclosing(values, value)


def test_run_reaction_calc():
    temperature = 311
    pressure = 2

    # Assuming this function returns a list with the Gibbs free energy values
    gibbs_values, equilibrium, reaction_ids = _interp_gibbs(temperature, pressure)

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


def test_interp_gibbs():
    # Mock data for reactions and pressure/temperature combinations
    reactions = [
        {1: {"e": None, "k": 1.0, "g": 100.0}, 2: {"e": None, "k": 1.5, "g": 150.0}},
        {1: {"e": None, "k": 1.0, "g": 100.0}, 2: {"e": None, "k": 1.5, "g": 150.0}},
        {1: {"e": None, "k": 1.0, "g": 100.0}, 2: {"e": None, "k": 1.5, "g": 150.0}},
        {1: {"e": None, "k": 2.0, "g": 200.0}, 2: {"e": None, "k": 2.5, "g": 250.0}},
    ]
    pt_combinations = [(350, 30), (350, 35), (400, 30), (400, 35)]
    pressure = 31
    temperature = 375

    # Call the function
    calculated_gibbs_values, equilibrium, reaction_ids = _interp_gibbs(
        reactions, pt_combinations, pressure, temperature
    )

    # Check if the returned values are of expected types
    assert isinstance(calculated_gibbs_values, list)
    assert isinstance(equilibrium, list)
    assert isinstance(reaction_ids, list)

    # Check if calculated_values matches expected length
    assert len(calculated_gibbs_values) == len(reaction_ids)

    # Check if equilibrium has the expected number of entries
    assert len(equilibrium) > 0


def test__interpolated_result_accuracy():
    for temp in range(len(TEMPERATURES) - 1):
        for press in range(len(PRESSURE_LIST) - 1):
            _interpolated_result_accuracy(temp, press)


def _interpolated_result_accuracy(temps, press):
    # Set desired temperature and pressure

    min_temp_idx = max(temps - 1, 0)
    max_temp_idx = min(temps + 1, len(TEMPERATURE_LIST) - 1)

    min_pressure_idx = max(press - 1, 0)
    max_pressure_idx = min(press + 1, len(PRESSURE_LIST) - 1)

    temperature = TEMPERATURE_LIST[temps]
    pressure = PRESSURE_LIST[press]

    # Get the true values to compare the results at the end
    true_reactions = get_reactions(temperature, pressure)

    # Collection reactions for enclosing Temperature and Pressure pairs
    # and sort the reactions in ascending order wrt reaction_ids

    # (T0, P0), (T0, P1), (T1, P0), (T1, P1)
    tp_combinations = [
        (TEMPERATURE_LIST[min_temp_idx], PRESSURE_LIST[min_pressure_idx]),
        (TEMPERATURE_LIST[min_temp_idx], PRESSURE_LIST[max_pressure_idx]),
        (TEMPERATURE_LIST[max_temp_idx], PRESSURE_LIST[min_pressure_idx]),
        (TEMPERATURE_LIST[max_temp_idx], PRESSURE_LIST[max_pressure_idx]),
    ]

    reactions: List[Dict[int, ReactionType]] = [
        get_reactions(t, p) for (t, p) in tp_combinations
    ]

    sorted_reactions = [dict(sorted(reaction.items())) for reaction in reactions]

    interpolated_gibbs_constant, _, ids = _interp_gibbs(
        sorted_reactions, tp_combinations, pressure, temperature
    )

    # Assertions to verify the accuracy of the interpolated results
    for reaction_id, interpolated_value in zip(ids, interpolated_gibbs_constant):
        true_value = true_reactions[reaction_id]["g"]
        assert np.isclose(interpolated_value, true_value, atol=0.5, rtol=0.01), (
            f"Interpolated value for reaction ID {reaction_id} "
            f"({interpolated_value}) does not match true value "
            f"({true_value})."
        )


def test_interpolate_gibbs_values():
    values = np.array([1269, 1519])
    point1 = 250
    point2 = 300
    target = 275
    expected_result = [1394]

    result = interpolate_gibbs_values(values, point1, point2, target)

    assert (expected_result == result).all()

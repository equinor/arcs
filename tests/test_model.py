import numpy as np
import pytest

from arcs.model import (
    _find_enclosing,
    get_reactions,
    _calculate_k,
    get_interpolated_reactions,
)

PRESSURE_LIST = [
    1,
    2,
    5,
    10,
    15,
    20,
    25,
    30,
    35,
    40,
    45,
    50,
    100,
    125,
    150,
    175,
    200,
    225,
    250,
    275,
    300,
]
TEMPERATURE_LIST = [200, 250, 300, 350, 400]


@pytest.mark.parametrize(
    "target, expected_lower, expected_upper",
    [
        (223.6, 200, 250),
        (200.0, 200, 250),
        (275.0, 250, 300),
        (399.9, 350, 400),
    ],
)
def test_find_enclosing(target: float, expected_lower: int, expected_upper: int):
    lower, upper = _find_enclosing(target, TEMPERATURE_LIST)
    assert lower == expected_lower
    assert upper == expected_upper


def test_interpolated_result_accuracy():
    for temp in range(len(TEMPERATURE_LIST) - 1):
        for press in range(len(PRESSURE_LIST) - 1):
            temperature = TEMPERATURE_LIST[temp]
            pressure = PRESSURE_LIST[press]

            true_reactions = get_reactions(temperature, pressure)
            interpolated_reactions = get_interpolated_reactions(temperature, pressure)

            assert [
                np.isclose(
                    interpolated_reactions[reaction_id]["g"],
                    true_reactions[reaction_id]["g"],
                )
                for reaction_id in range(len(true_reactions))
            ]


def test_calculate_k():
    for temperature in TEMPERATURE_LIST:
        for pressure in PRESSURE_LIST:
            reactions_list = get_reactions(temperature, pressure)
            for reaction_id in range(len(reactions_list)):
                gibbs_energy = reactions_list[reaction_id]["g"]
                true_k = reactions_list[reaction_id]["k"]
                calculated_k = _calculate_k(gibbs_energy, temperature)
                assert np.isclose(calculated_k, true_k, rtol=1e-4, atol=1e-4)

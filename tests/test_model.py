import numpy as np
import pytest
import copy

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
        (200.0, 200, 200),
        (275.0, 250, 300),
        (399.9, 350, 400),
    ],
)
def test_find_enclosing(target: float, expected_lower: int, expected_upper: int):
    enclosing = _find_enclosing(target, TEMPERATURE_LIST)
    assert enclosing[0] == expected_lower
    assert enclosing[-1] == expected_upper


def test_interpolated_handle_points_in_list():
    # We only assert the structure of the return type
    result = get_interpolated_reactions(TEMPERATURE_LIST[1], PRESSURE_LIST[1])

    assert len(result) > 9000
    assert isinstance(result, dict)
    assert all(x in result[0] for x in ["e", "g", "k"])


def test_interpolated_handle_single_axis():
    # interpolation along temperature
    for temperature in TEMPERATURE_LIST[1:-1]:
        for pressure in PRESSURE_LIST[0:-1]:
            mocked_temperatures = copy.copy(TEMPERATURE_LIST)
            mocked_temperatures.remove(temperature)

            with pytest.MonkeyPatch.context() as m:
                m.setattr("arcs.model.TEMPERATURE_LIST", mocked_temperatures)

                true_reactions = get_reactions(temperature, pressure)
                interpolated_reactions = get_interpolated_reactions(
                    temperature, pressure
                )

                assert [
                    np.isclose(
                        interpolated_reactions[reaction_id]["g"],
                        true_reactions[reaction_id]["g"],
                    )
                    for reaction_id in range(len(true_reactions))
                ]

    # interpolation along pressure
    for temperature in TEMPERATURE_LIST[0:-1]:
        for pressure in PRESSURE_LIST[1:-1]:
            mocked_pressures = copy.copy(PRESSURE_LIST)
            mocked_pressures.remove(pressure)

            with pytest.MonkeyPatch.context() as m:
                m.setattr("arcs.model.PRESSURE_LIST", mocked_pressures)

                true_reactions = get_reactions(temperature, pressure)
                interpolated_reactions = get_interpolated_reactions(
                    temperature, pressure
                )

                assert [
                    np.isclose(
                        interpolated_reactions[reaction_id]["g"],
                        true_reactions[reaction_id]["g"],
                    )
                    for reaction_id in range(len(true_reactions))
                ]


def test_interpolated_result_accuracy():
    for temperature in TEMPERATURE_LIST[1:-1]:
        for pressure in PRESSURE_LIST[1:-1]:
            mocked_temperatures = copy.copy(TEMPERATURE_LIST)
            mocked_pressures = copy.copy(PRESSURE_LIST)
            mocked_temperatures.remove(temperature)
            mocked_pressures.remove(pressure)

            with pytest.MonkeyPatch.context() as m:
                m.setattr("arcs.model.PRESSURE_LIST", mocked_pressures)
                m.setattr("arcs.model.TEMPERATURE_LIST", mocked_temperatures)

                true_reactions = get_reactions(temperature, pressure)
                interpolated_reactions = get_interpolated_reactions(
                    temperature, pressure
                )

                assert [
                    np.isclose(
                        interpolated_reactions[reaction_id]["g"],
                        true_reactions[reaction_id]["g"],
                    )
                    for reaction_id in range(len(true_reactions))
                ]


def test_calculate_k():
    get_reactions.cache_clear()
    for temperature in TEMPERATURE_LIST:
        for pressure in PRESSURE_LIST:
            reactions_list = get_reactions(temperature, pressure)
            for reaction_id in range(len(reactions_list)):
                gibbs_energy = reactions_list[reaction_id]["g"]
                true_k = reactions_list[reaction_id]["k"]
                calculated_k = _calculate_k(gibbs_energy, temperature)
                assert np.isclose(calculated_k, true_k)

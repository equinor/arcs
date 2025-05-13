from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import pickle
from typing import List, Tuple, TypedDict, Dict, Any
import chempy
import pandas as pd
import numpy as np
import numpy.typing as npt
from scipy import interpolate
from scipy.interpolate import LinearNDInterpolator

from bin.generate_tables import process_generic_inputs

MODEL_PATH = Path(__file__).parent.parent / "model"

BOLTZMANN_CONSTANT = 8.617333262 * (10 ** (-5))  # Boltzmann constant k in eV

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


class ReactionType(TypedDict):
    e: chempy.Equilibrium
    k: float
    g: float


@lru_cache
def get_reactions(temperature: int, pressure: int) -> dict[int, ReactionType]:
    with open(
        MODEL_PATH / f"T{temperature}_P{pressure}" / "reactions.p", "rb"
    ) as stream:
        return pickle.load(stream)  # type: ignore


Table = dict[tuple[str, str], tuple[npt.NDArray[np.int16], npt.NDArray[np.float64]]]


@lru_cache
def get_table(
    temperature: int, pressure: int, *, rank_small_reactions_higher: bool = True
) -> Table:
    suffix = "-rank_small_reactions_higher" if rank_small_reactions_higher else ""
    file_path = MODEL_PATH / f"T{temperature}_P{pressure}" / f"table{suffix}.p"

    if not file_path.exists():
        reactions = get_interpolated_reactions(temperature, pressure)
        process_generic_inputs(
            reactions, temperature, pressure, MODEL_PATH
        )
        destination_directory = MODEL_PATH / f"T{temperature}_P{pressure}"
        reactions_file_path = destination_directory / "reactions.p"
        with open(reactions_file_path, "wb") as stream:
            pickle.dump(reactions, stream)

    with open(file_path, "rb") as stream:
        return pickle.load(stream)


def get_interpolated_reactions(temperature: float, pressure: float) -> dict[int, ReactionType]:
    pressure_list: Tuple[int, int] = _find_enclosing(pressure, PRESSURE_LIST)
    temperature_list: Tuple[int, int] = _find_enclosing(temperature, TEMPERATURE_LIST)

    results = []
    tlogp_combinations = [[t, np.log(p)] for t in temperature_list for p in pressure_list]
    for p in pressure_list:
        for t in temperature_list:
            reactions = get_reactions(t, p)
            for reaction_id in range(len(reactions)):
                g_value = reactions[reaction_id]['g']
                results.append({
                    f"Reaction_id": reaction_id,
                    f"Gibbs_{t}_{p}": g_value
                })
    df_results = pd.DataFrame(results)
    df_grouped = df_results.groupby("Reaction_id").agg(
        lambda x: x.dropna().iloc[0] if not x.dropna().empty else None).reset_index()

    interpolators = [LinearNDInterpolator(tlogp_combinations, np.array(df_grouped.iloc[i, 1:df_grouped.columns.size])) for i in range(len(df_grouped))]
    gibbs_values = [interpolator(temperature, np.log(pressure)) for interpolator in interpolators]
    reactions_table = {}
    for reaction_id in range(len(reactions)):
        reaction: ReactionType = reactions[reaction_id]
        reaction['k'] = _calculate_k(gibbs_values[reaction_id], temperature)
        reaction['g'] = gibbs_values[reaction_id]
        reactions_table[reaction_id] = reaction

    return reactions_table


def _calculate_k(gibbs_energy: float, temperature: float) -> np.float64:
    try:
        return np.float64(
            np.exp((-1) * gibbs_energy / (BOLTZMANN_CONSTANT * temperature))
        )
    except (OverflowError, ValueError):
        return np.float64(0)


def run_reaction_calc(
    pressure: int, temperature: int
) -> Tuple[List[float], List[chempy.Equilibrium], List[int]]:
    pressure_list = [
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
    temperature_list = [200, 250, 300, 350, 400]

    temperature_boundaries: Tuple[int, int] = _find_enclosing(
        temperature, temperature_list
    )
    pressure_boundaries: Tuple[int, int] = _find_enclosing(pressure, pressure_list)

    pt_combinations: List[Tuple[int, int]] = [
        (t, p) for t in temperature_boundaries for p in pressure_boundaries
    ]

    reactions: List[Dict[int, ReactionType]] = [
        get_reactions(t, p) for (t, p) in pt_combinations
    ]

    sorted_reactions = [dict(sorted(reaction.items())) for reaction in reactions]

    return _get_gibbs_constant(sorted_reactions, pt_combinations, pressure, temperature)


def _get_gibbs_constant(
    enclosing_reactions: List[Dict[int, ReactionType]],
    pt_combinations: List[Tuple[int, int]],
    pressure: int,
    temperature: int,
) -> Tuple[List[float], List[chempy.Equilibrium], List[int]]:
    """
    Calculates Gibbs contant using linear interpolation between 2 enclosing points

    Returns: Data required for reaction table
    """

    gibbs_values: List[float] = []
    equilibrium: List[chempy.Equilibrium] = []

    for reactions in enclosing_reactions:  # Runs for each P & T enclosing set
        reaction_gibbs_values = []
        reaction_ids = []
        for (
            reaction_id,
            reaction_data,
        ) in reactions.items():  # Runs for all 9113 reaction ids
            g_value = reaction_data["g"]
            _e = reaction_data["e"]
            equilibrium.append(_e)
            reaction_gibbs_values.append(g_value)
            reaction_ids.append(reaction_id)
        gibbs_values.append(reaction_gibbs_values)  # type: ignore

    gibbs_values = np.array(gibbs_values)  # type: ignore

    # Unpack list of tuples
    (T0, P0), (T0, P1), (T1, P0), (T1, P1) = pt_combinations

    # Keeping pressure constant
    gibbs_values_p_low = np.array([gibbs_values[0], gibbs_values[2]]).transpose()
    gibbs_values_p_high = np.array([gibbs_values[1], gibbs_values[3]]).transpose()

    # Interpolate Gibbs values for temperature
    g_val_low_pressure = interpolate_gibbs_values(
        gibbs_values_p_low, T0, T1, temperature
    )
    g_val_high_pressure = interpolate_gibbs_values(
        gibbs_values_p_high, T0, T1, temperature
    )

    # Prepare for pressure interpolation
    gibbs_values_for_changing_pressure = np.array(
        [g_val_low_pressure, g_val_high_pressure]
    ).transpose()

    # Interpolate Gibbs values for pressure
    calculated_gibbs_values = interpolate_gibbs_values(
        gibbs_values_for_changing_pressure, P0, P1, pressure
    )

    return [calculated_gibbs_values, equilibrium, reaction_ids]  # type: ignore


def interpolate_gibbs_values(
    values: npt.NDArray[np.float64], point1: int, point2: int, target: int
) -> List[float]:
    interpolation = [interpolate.interp1d([point1, point2], y) for y in values]
    return [f(target) for f in interpolation]


def _find_enclosing(target: float, values: List[int]) -> Tuple[int, int]:
    lower_boundary = max(max([x for x in values if x <= target]), values[0])
    upper_boundary = min(min([x for x in values if x > target]), values[-1])
    return lower_boundary, upper_boundary


def get_reaction_compounds(reactions: dict[int, ReactionType]) -> dict[int, set[str]]:
    return {k: set(r["e"].reac) | set(r["e"].prod) for k, r in reactions.items()}

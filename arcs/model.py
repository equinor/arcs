from __future__ import annotations

import pickle
from functools import lru_cache
from pathlib import Path
from typing import List, TypedDict

import chempy
import numpy as np
import numpy.typing as npt
from scipy.interpolate import LinearNDInterpolator  # type: ignore

from bin.generate_tables import process_generic_inputs

MODEL_PATH = Path(__file__).parent.parent / "model"

BOLTZMANN_CONSTANT = 8.617333262 * (10 ** (-5))  # Boltzmann constant k in eV

TEMPERATURE_LIST = [200, 250, 300, 350, 400]
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
        process_generic_inputs(reactions, temperature, pressure, MODEL_PATH)
        destination_directory = MODEL_PATH / f"T{temperature}_P{pressure}"
        reactions_file_path = destination_directory / "reactions.p"
        with open(reactions_file_path, "wb") as stream:
            pickle.dump(reactions, stream)

    with open(file_path, "rb") as stream:
        return pickle.load(stream)  # type: ignore


def get_interpolated_reactions(
    temperature: float, pressure: float
) -> dict[int, ReactionType]:
    temperature_boundary = _find_enclosing(temperature, TEMPERATURE_LIST)
    pressure_boundary = _find_enclosing(pressure, PRESSURE_LIST)

    const_temperature = temperature_boundary[0] == temperature_boundary[-1]
    const_pressure = pressure_boundary[0] == pressure_boundary[-1]

    if const_temperature and const_pressure:
        # no need to do interpolation
        return get_reactions(temperature, pressure)

    results: dict[int, List[float]] = {}
    tlogp_combinations = [
        [t, np.log(p)] for t in temperature_boundary for p in pressure_boundary
    ]
    for t in temperature_boundary:
        for p in pressure_boundary:
            reactions = get_reactions(t, p)
            for reaction_id in range(len(reactions)):
                g_value = reactions[reaction_id]["g"]
                if reaction_id not in results:
                    results[reaction_id] = []
                results[reaction_id].append(g_value)

    if const_pressure:
        # interpolate along temperature axis
        t_vals = [t[0] for t in tlogp_combinations]
        gibbs_values = [
            np.interp(temperature, t_vals, results[i]) for i in range(len(results))
        ]

    if const_temperature:
        # interpolate along pressure axis
        p_vals = [p[1] for p in tlogp_combinations]
        gibbs_values = [
            np.interp(pressure, p_vals, results[i]) for i in range(len(results))
        ]

    if not const_temperature and not const_pressure:
        # interpolate both axis
        interpolators = [
            LinearNDInterpolator(tlogp_combinations, np.array(results[i]))
            for i in range(len(results))
        ]

        gibbs_values = [
            interpolator(temperature, np.log(pressure))
            for interpolator in interpolators
        ]

    reactions_table = {}
    for reaction_id in range(len(reactions)):
        reaction: ReactionType = reactions[reaction_id]
        k = _calculate_k(gibbs_values[reaction_id], temperature)
        g = gibbs_values[reaction_id]
        reactions_table[reaction_id] = ReactionType(e=reaction["e"], k=k, g=g)

    return reactions_table


def _calculate_k(gibbs_energy: float, temperature: float) -> np.float64:
    try:
        return np.float64(
            np.exp((-1) * gibbs_energy / (BOLTZMANN_CONSTANT * temperature))
        )
    except (OverflowError, ValueError):
        return np.float64(0)


def _find_enclosing(target: float, values: List[int]) -> tuple[int, int]:
    lower_boundary = max(max([x for x in values if x <= target]), values[0])
    upper_boundary = min(min([x for x in values if x >= target]), values[-1])
    return [lower_boundary, upper_boundary]


def get_reaction_compounds(reactions: dict[int, ReactionType]) -> dict[int, set[str]]:
    return {k: set(r["e"].reac) | set(r["e"].prod) for k, r in reactions.items()}

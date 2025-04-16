from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import pickle
import shutil
from typing import List, Tuple, TypedDict, Dict
import chempy
import numpy as np
import numpy.typing as npt
from scipy import interpolate
from bin.generate_tables import process_generic_inputs

MODEL_PATH = Path(__file__).parent.parent / "model"


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
        
        g,e,id = run_reaction_calc(pressure, temperature)
        restructured_reactions = [
            {_id:{
                "e":_e,
                "k": np.nan,# Oneline of FunctionCall(),
                "g":_g ,        
            }}
            for _id, _e, _g in zip(id, e, g)]
        
        process_generic_inputs(restructured_reactions, temperature, pressure, MODEL_PATH)
        
        source_file = MODEL_PATH / "T200_P1" / "reactions.p"
        destination_directory = MODEL_PATH / f"T{temperature}_P{pressure}"
        shutil.copy(source_file, destination_directory / "reactions.p")

    with open(file_path, "rb") as stream:
        return pickle.load(stream)  # type: ignore


def run_reaction_calc(pressure: int, temperature: int):
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

    return _get_gibbs_constant(reactions, pt_combinations, pressure, temperature)


def _get_gibbs_constant(
    reactions: List[Dict[int, ReactionType]],
    pt_combinations: List[Tuple[int, int]],
    pressure: int,
    temperature: int,
) -> List[np.float64]:
    """
    Input Arguments
    pt_combinations = [(T0,P0), (T0, P1), (T1,P0), (T1,P1)]
    reactions = [{reaction_id_1: {e : chempy_obj ,k : np.float64 ,g: np.float64}, reaction_id_1: {e,k,g} }, ... ]

    """

    gibbs_values = []
    equilibrium = []

    for reaction in reactions: # - runs 4 times because we have 4 combinations of P and T
        reaction_gibbs_values = []
        reaction_ids = []
        for reaction_id, reaction_data in reaction.items(): # runs 9113 times 
            g_value = reaction_data["g"]
            _e = reaction_data["e"]
            equilibrium.append(_e) 
            reaction_gibbs_values.append(g_value)
            reaction_ids.append(reaction_id)

        gibbs_values.append(reaction_gibbs_values)

    gibbs_values = np.array(gibbs_values) #[4 numpy arrays - with 9113 elements 4 x 9113, but i need transpose of this]

    # Unpack list of tuples
    (T0, P0), (T0, P1), (T1, P0), (T1, P1) = pt_combinations

    # Keeping pressure constant
    gibbs_values_p_low = np.array([gibbs_values[0], gibbs_values[2]]).transpose()
    gibbs_values_p_high = np.array([gibbs_values[1], gibbs_values[3]]).transpose()

    # Linear interpolate for temeprature
    interpolation_g_low_pressure = [
        interpolate.interp1d([T0, T1], y) for y in gibbs_values_p_low
    ]
    interpolation_g_high_pressure = [
        interpolate.interp1d([T0, T1], y) for y in gibbs_values_p_high
    ]

    g_val_low_pressure = [g(temperature) for g in interpolation_g_low_pressure]
    g_val_high_pressure = [g(temperature) for g in interpolation_g_high_pressure]

    # Logarithmic interpolation on pressures
    gibbs_values_for_changing_pressure = np.array(
        [g_val_low_pressure, g_val_high_pressure]
    ).transpose()
    interpolation_g = [
        interpolate.interp1d([P0, P1], y) for y in gibbs_values_for_changing_pressure
    ]

    calculated_gibbs_values = [g(pressure) for g in interpolation_g]

    return [calculated_gibbs_values, equilibrium, reaction_ids]


def _find_enclosing(target: int, values: List[int]):

    if (target in values):
        return target, target
    lower_candidate = max(list(filter(lambda x: x < target, values)))
    upper_candidate = min(list(filter(lambda x: x > target, values)))

    return lower_candidate, upper_candidate


def get_reaction_compounds(reactions: dict[int, ReactionType]) -> dict[int, set[str]]:
    return {k: set(r["e"].reac) | set(r["e"].prod) for k, r in reactions.items()}

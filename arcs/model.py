from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import pickle
from typing import TypedDict
import chempy
import numpy as np
import numpy.typing as npt


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

    with open(
        MODEL_PATH / f"T{temperature}_P{pressure}" / f"table{suffix}.p", "rb"
    ) as stream:
        return pickle.load(stream)  # type: ignore


def get_reaction_compounds(reactions: dict[int, ReactionType]) -> dict[int, set[str]]:
    return {k: set(r["e"].reac) | set(r["e"].prod) for k, r in reactions.items()}

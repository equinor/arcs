from __future__ import annotations
from functools import lru_cache
from pathlib import Path
import pickle
from typing import TypedDict
import chempy
import networkx as nx


MODEL_PATH = Path(__file__).parent.parent / "model"


class ReactionType(TypedDict):
    e: chempy.Equilibrium
    k: float
    g: float


@lru_cache
def get_graph(temperature: int, pressure: int) -> nx.MultiDiGraph:
    with open(MODEL_PATH / f"T{temperature}_P{pressure}" / "graph.p", "rb") as stream:
        return pickle.load(stream)


@lru_cache
def get_reactions(temperature: int, pressure: int) -> dict[int, ReactionType]:
    with open(
        MODEL_PATH / f"T{temperature}_P{pressure}" / "reactions.p", "rb"
    ) as stream:
        return pickle.load(stream)

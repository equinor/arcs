#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
from collections import defaultdict
import os
from typing import IO, TYPE_CHECKING, Any
import numpy as np
import pickle
import re
from multiprocessing import Pool


if TYPE_CHECKING:
    from arcs.model import ReactionType, Table


MOLECULE_PATTERN = re.compile(r"[A-Z][a-z]?(\d+)?")


def _num_atoms(name: str) -> int:
    return sum(int(x) if x else 1 for x in MOLECULE_PATTERN.findall(name))


def _cost_function(
    gibbs: float,
    temperature: float,
    compounds: dict[str, int],
) -> float:
    num_atoms = sum(_num_atoms(k) * v for k, v in compounds.items())
    return np.log(1 + (273 / temperature) * np.exp(gibbs / num_atoms))


def _compute_edge(
    start: str,
    gibbs: float,
    temperature: float,
    reac: dict[str, int],
    prod: dict[str, int],
) -> float:
    assert (
        start in reac or start in prod
    ), f"{start} is neither in reactants ({reac}) or products ({prod})"

    if start in reac:
        return _cost_function(gibbs, temperature, reac)
    elif start in prod:
        return _cost_function(-gibbs, temperature, prod)


_PreTable = dict[tuple[str, str], list[tuple[int, float, float]]]


def _process_reaction(
    table: _PreTable, temperature: float, index: int, reaction: ReactionType
) -> None:
    eq = reaction["e"]
    components: set[str] = set(eq.reac) | set(eq.prod)

    for src in components:
        for dst in components:
            if src == dst:
                continue

            cost = _compute_edge(src, reaction["g"], temperature, eq.reac, eq.prod)

            # For when `rank_small_reactions_higher is True`
            cost_small = cost * 10 ** len(components)

            table[(src, dst)].append((index, cost, cost_small))


def _write_table(
    stream: IO[bytes], pre_table: _PreTable, value_column: int = 1
) -> None:
    table: Table = {}
    for path, data in pre_table.items():
        reaction_ids = np.array([x[0] for x in data], dtype=np.int16)
        costs = np.array([x[value_column] for x in data], dtype=np.float64)

        # Sort costs in decreasing order
        indices = np.argsort(costs, kind="stable")
        reaction_ids = reaction_ids[indices]
        costs = 1 / costs[indices]

        table[path] = (reaction_ids, costs)
    pickle.dump(table, stream)


def _process_file(filepath: Path) -> None:
    print(f"Processing: {filepath}   ", end="\r")
    with filepath.open("rb") as f:
        reactions: dict[int, ReactionType] = pickle.load(f)

    match = re.match(r"T(\d+)_P\d+", filepath.parent.name)
    assert (
        match is not None
    ), f"Couldn't get temperature information from '{filepath.parent.name}' model"

    temperature = int(match[1])

    pre_table: _PreTable = defaultdict(list)
    for index, reaction in reactions.items():
        _process_reaction(pre_table, temperature, index, reaction)

    with (filepath.parent / "table.p").open("wb") as f:
        _write_table(f, pre_table)
    with (filepath.parent / "table-rank_small_reactions_higher.p").open("wb") as f:
        _write_table(f, pre_table, value_column=2)


def process_generic_inputs(
    reactions: dict[Any, dict[str, Any]],
    temperature: int,
    pressure: int,
    path: Path,
):
    directory_name = path / f"T{temperature}_P{pressure}"

    try:
        os.makedirs(directory_name, exist_ok=True)
        print("Directory created successfully", directory_name)

    except Exception as e:
        print(f"An error occurred: {e}")

    pre_table: _PreTable = defaultdict(list)

    for index, values in reactions.items():
        _process_reaction(pre_table, temperature, index, values)

    with (directory_name / "table.p").open("wb") as f:
        _write_table(f, pre_table)

    with (directory_name / "table-rank_small_reactions_higher.p").open("wb") as f:
        _write_table(f, pre_table, value_column=2)


def main() -> None:
    with Pool() as pool:
        _ = pool.map(_process_file, Path().glob("model/*/reactions.p"))
    print("\nDone!")


if __name__ == "__main__":
    main()

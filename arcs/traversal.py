from __future__ import annotations
from collections import defaultdict
from collections.abc import Collection
from dataclasses import dataclass
from typing import Literal, Any, TypedDict

from chempy.equilibria import Equilibrium, EqSystem
from chempy import Substance
import copy
import networkx as nx
import itertools as it
import concurrent.futures
import platform
import psutil

from datetime import datetime
import numpy as np
import pandas as pd

from arcs.model import get_graph, get_reactions
import warnings

warnings.filterwarnings("ignore")


@dataclass
class TraversalResult:
    initfinaldiff: Any
    final_concs: Any
    metadata: dict[str, Any]
    data: dict[int, Any]

    def to_dict(self) -> dict[Any, Any]:
        return {
            "initfinaldiff": self.initfinaldiff,
            "final_concs": self.final_concs,
            "metadata": self.metadata,
            "data": {k: v for k, v in self.data.items()},
        }


def _get_weighted_random_compounds(
    temperature: int,
    pressure: int,
    init_concs: dict[str, float],
    *,
    co2: bool,
    max_compounds: int,
    probability_threshold: float,
    scale_highest: float,
    ceiling: int,
    rng: np.random.Generator,
) -> dict[Any, Any]:
    concs = copy.deepcopy(init_concs)
    if co2 is False and "CO2" in concs:
        del concs["CO2"]

    num_positive_concs = len([conc for conc in concs.values() if conc > 0])
    if max_compounds > num_positive_concs:
        max_compounds = num_positive_concs

    median_conc = np.median([conc for conc in concs.values() if conc > 0])
    concs_above_ceiling = {
        k: v for k, v in concs.items() if v > (median_conc * (1 + (ceiling / 100)))
    }

    for compound, conc in concs_above_ceiling.items():
        concs[compound] = conc * scale_highest

    compound_probabilities = {k: v / sum(concs.values()) for k, v in concs.items()}
    filtered_probabilities = {
        k: v for k, v in compound_probabilities.items() if v > probability_threshold
    }
    compound_probabilities = {
        k: v / sum(filtered_probabilities.values())
        for k, v in filtered_probabilities.items()
    }

    sample_frame = list(
        rng.choice(
            list(compound_probabilities.keys()),
            100,
            p=list(compound_probabilities.values()),
        )
    )
    selected_compounds = {}
    for c in range(max_compounds):
        if c == 0:
            choice = rng.choice(sample_frame)
            selected_compounds[choice] = compound_probabilities[choice]
        else:
            try:
                for i in range(sample_frame.count(list(selected_compounds)[c - 1])):
                    sample_frame.remove(list(selected_compounds)[c - 1])
                try:
                    choice = rng.choice(sample_frame)
                    selected_compounds[choice] = compound_probabilities[choice]
                except Exception:
                    pass
            except Exception:
                pass

    return selected_compounds


def _length_multiplier(
    candidate: list[Any], *, rank_small_reactions_higher: bool
) -> int:
    if rank_small_reactions_higher:
        return len(list(candidate))
    else:
        return 1


def _get_weighted_reaction_rankings(
    tempreature: int,
    pressure: int,
    choices: Collection[str],
    *,
    max_rank: int,
    method: Literal["Bellman-Ford", "Dijkstra"],
    rank_small_reactions_higher: bool,
    graph: nx.MultiDiGraph,
) -> dict[Any, dict[str, Any]] | None:
    reaction_rankings = {}
    if len(choices) > 1:
        possible_shortest_reaction_paths = list(
            nx.shortest_paths.all_shortest_paths(
                graph, list(choices)[0], list(choices)[1], method=method
            )
        )
        for reaction_path in possible_shortest_reaction_paths:
            reactant_compound = reaction_path[0]
            reaction = reaction_path[1]
            candidates = list(graph[reaction])
            if len(choices) > 2:
                for compound in list(choices)[2:]:
                    if compound in candidates:
                        weight = graph.get_edge_data(reactant_compound, reaction)[0][
                            "weight"
                        ] * 10 ** _length_multiplier(
                            graph[reaction],
                            rank_small_reactions_higher=rank_small_reactions_higher,
                        )
                        reaction_rankings[reaction] = {
                            "candidates": candidates,
                            "weight": weight,
                        }
            else:
                weight = graph.get_edge_data(reactant_compound, reaction)[0][
                    "weight"
                ] * 10 ** _length_multiplier(
                    graph[reaction],
                    rank_small_reactions_higher=rank_small_reactions_higher,
                )
                reaction_rankings[reaction] = {
                    "candidates": candidates,
                    "weight": weight,
                }
    if reaction_rankings:
        sorted_rankings = (
            pd.DataFrame(reaction_rankings).sort_values(by="weight", axis=1).to_dict()
        )
        top_ranked_reactions = [
            reaction for i, reaction in enumerate(sorted_rankings) if i <= max_rank
        ]
        reaction_rankings = {r: reaction_rankings[r] for r in top_ranked_reactions}
        return reaction_rankings
    else:
        return None


def _generate_eqsystem(
    index: int, temperature: int, pressure: int, *, reactions: dict[int, Any]
) -> EqSystem | None:
    charged_species = {
        "CO3H": -1,
        "NH4": +1,
        "NH2CO2": -1,
    }
    reaction = reactions[index]
    reactants = reaction["e"].reac
    products = reaction["e"].prod
    k = reaction["k"]
    substances = {}

    reaction_compounds = list(it.chain(*[list(reactants) + list(products)]))
    for compound in reaction_compounds:
        if compound in list(charged_species.keys()):
            substance = Substance.from_formula(
                compound, **{"charge": charged_species[compound]}
            )
            substances[substance.name] = substance
        else:
            substance = Substance.from_formula(
                compound, **{"charge": 0}
            )  # ,charge=0) # charge buggers up everything have removed for now....
            substances[substance.name] = substance
    equilibrium = Equilibrium(reac=reactants, prod=products, param=k)
    try:
        return EqSystem(
            [equilibrium], substances
        )  # might not just be able to try a return...
    except Exception:
        return None


def _equilibrium_concentrations(
    concs: dict[str, float], eq: EqSystem
) -> tuple[dict[str, float], str]:
    # something is going wrong here...
    equilibrium_concentrations = defaultdict(lambda: 0.0, concs)
    try:
        root_concs, solution_data, sane = eq.root(equilibrium_concentrations)
        assert solution_data["success"] and sane
        for i, conc in enumerate(root_concs):
            equilibrium_concentrations[eq.substance_names()[i]] = conc
        concs = equilibrium_concentrations
        eq = eq.string()
    except Exception:
        concs = equilibrium_concentrations
        eq = None
    return (dict(concs), eq)


class _RandomWalk(TypedDict):
    data: dict[str, float]
    equation_statistics: list[Equilibrium]
    path_length: int


def _random_walk(
    temperature: int,
    pressure: int,
    concs: dict[str, float],
    *,
    probability_threshold: float,
    path_depth: int,
    max_compounds: int,
    max_rank: int,
    co2: bool,
    scale_highest: float,
    ceiling: int,
    method: Literal["Bellman-Ford", "Dijkstra"],
    rank_small_reactions_higher: bool,
    rng: np.random.Generator,
    reactions: dict[int, Any],
    graph: nx.MultiDiGraph,
) -> _RandomWalk:
    conc_history = [concs]
    reaction_history: list[str] = []

    for _ in range(path_depth):
        previous_conc_step = conc_history[-1]
        try:
            choices = _get_weighted_random_compounds(
                temperature,
                pressure,
                previous_conc_step,
                max_compounds=max_compounds,
                probability_threshold=probability_threshold,
                co2=co2,
                scale_highest=scale_highest,
                ceiling=ceiling,
                rng=rng,
            )
        except Exception:
            break
        if len(choices) <= 1:
            break
        reaction_rankings = _get_weighted_reaction_rankings(
            temperature,
            pressure,
            choices,
            max_rank=max_rank,
            method=method,
            rank_small_reactions_higher=rank_small_reactions_higher,
            graph=graph,
        )
        if not reaction_rankings:
            break
        reaction_weights = {
            r: 1 / reaction_rankings[r]["weight"] for r in reaction_rankings
        }
        probabilities = {
            k: v / sum(reaction_weights.values()) for k, v in reaction_weights.items()
        }
        chosen_reaction = rng.choice(
            [
                rng.choice(
                    list(probabilities.keys()),
                    len(probabilities),
                    p=list(probabilities.values()),
                )
            ][0]
        )

        eqsyst = _generate_eqsystem(
            chosen_reaction, temperature, pressure, reactions=reactions
        )
        # if reaction was previous reaction then break
        path_available = [r for r in reaction_history if r is not None]
        if path_available:
            assert eqsyst is not None
            if (
                eqsyst.string() == path_available[-1]
                and eqsyst.string() == path_available[-1]
            ):
                break  # extra break

        equilibrium_concs, reaction_string = _equilibrium_concentrations(
            previous_conc_step, eqsyst
        )
        conc_history.append(equilibrium_concs)
        reaction_history.append(reaction_string)

    return {
        "data": conc_history[-1],
        "equation_statistics": [r for r in reaction_history if r is not None],
        "path_length": len([r for r in reaction_history if r is not None]),
    }


def _sample_chunk(
    chunk_id: int, chunk_length: int, config: dict[str, Any]
) -> dict[int, _RandomWalk]:
    result_dict: dict[int, _RandomWalk] = {}
    for sample in range(chunk_length):
        result_dict[1 + sample + chunk_id * chunk_length] = _random_walk(**config)
        # +1 because we don't want to overwrite initial concs
    return result_dict


def _sample(
    temperature: int,
    pressure: int,
    concs: dict[str, float],
    *,
    co2: bool,
    max_compounds: int,
    probability_threshold: float,
    max_rank: int,
    sample_length: int,
    path_depth: int,
    ceiling: int,
    scale_highest: float,
    rank_small_reactions_higher: bool,
    method: Literal["Bellman-Ford", "Dijkstra"],
    rng: np.random.Generator,
    reactions: dict[int, Any],
    graph: nx.MultiDiGraph,
    nproc: int,
) -> dict[int, Any]:
    config = {
        "temperature": temperature,
        "pressure": pressure,
        "concs": concs,
        "co2": co2,
        "max_compounds": max_compounds,
        "probability_threshold": probability_threshold,
        "max_rank": max_rank,
        "path_depth": path_depth,
        "ceiling": ceiling,
        "scale_highest": scale_highest,
        "rank_small_reactions_higher": rank_small_reactions_higher,
        "method": method,
        "rng": rng,
        "reactions": reactions,
        "graph": graph,
    }
    if nproc == 0:
        nproc = psutil.cpu_count() or 1  # because cpu_count can return None

    chunk_size = sample_length // nproc

    result_dict: dict[int, _RandomWalk] = {
        0: {"data": concs, "equation_statistics": [], "path_length": 0}
    }

    if nproc == 1 or sample_length < nproc:
        print("Single process")
        result_dict.update(_sample_chunk(0, sample_length, config))
    else:
        print(f"CPU count: {nproc}, chunk size: {chunk_size}")
        with concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as executor:
            futures = {
                executor.submit(_sample_chunk, chunk_id, chunk_size, config): chunk_id
                for chunk_id in range(nproc)
            }
            for future in concurrent.futures.as_completed(futures):
                try:
                    result_dict.update(future.result())
                except Exception as exc:
                    chunk_id = futures[future]
                    print(f"Chunk {chunk_id} generated an exception: {exc}")

    return result_dict


def traverse(
    temperature: int,
    pressure: int,
    concs: dict[str, float],
    *,
    co2: bool = False,
    max_compounds: int = 5,
    probability_threshold: float = 0.05,
    max_rank: int = 5,
    sample_length: int = 1000,
    path_depth: int = 5,
    ceiling: int = 500,
    scale_highest: float = 0.1,
    rank_small_reactions_higher: bool = True,
    method: Literal["Bellman-Ford", "Dijkstra"] = "Dijkstra",
    rng: np.random.Generator | None = None,
    reactions: dict[int, Any] | None = None,
    graph: nx.MultiDiGraph | None = None,
    nproc: int = 0,
) -> TraversalResult:
    if rng is None:
        rng = np.random.default_rng()
    if graph is None:
        graph = get_graph(temperature, pressure)
    if reactions is None:
        reactions = get_reactions(temperature, pressure)

    concstring = pd.Series({k: v for k, v in concs.items() if v > 0}) / 1e-6
    if "CO2" in concstring:
        del concstring["CO2"]

    path_lengths = []
    samples = _sample(
        temperature,
        pressure,
        concs,
        co2=co2,
        max_compounds=max_compounds,
        probability_threshold=probability_threshold,
        max_rank=max_rank,
        sample_length=sample_length,
        path_depth=path_depth,
        ceiling=ceiling,
        scale_highest=scale_highest,
        rank_small_reactions_higher=rank_small_reactions_higher,
        method=method,
        rng=rng,
        reactions=reactions,
        graph=graph,
        nproc=nproc,
    )

    mean_concs = pd.DataFrame([s["data"] for s in samples.values()]).mean()
    df_summary = pd.DataFrame({"initial": concs, "final": mean_concs}) * 1e6
    df_summary = df_summary.dropna(how="all").fillna(0.0)
    df_summary["change"] = df_summary["final"] - df_summary["initial"]
    df_summary = df_summary.loc[(df_summary.abs() >= 1e-6).any(axis=1)]

    avg_path_length = np.median(
        [s["path_length"] for s in samples.values() if s["path_length"] is not None]
    )
    path_lengths.append(avg_path_length)

    metadata = {
        "arcs_version": "1.4.0",
        "avg_path_length": np.mean(path_lengths),
        "co2": co2,
        "max_compounds": max_compounds,
        "probability_threshold": probability_threshold,
        "shortest_path_method": method,
        "max_rank": max_rank,
        "sample_length": sample_length,
        "path_depth": path_depth,
        "ceiling": ceiling,
        "scale_highest": scale_highest,
        "rank_small_reactions_higher": rank_small_reactions_higher,
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "processor": platform.processor(),
        "available_cores": psutil.cpu_count(),
        "available_memory": str(int(psutil.virtual_memory()[0] / 1000 / 1000 / 1000))
        + "Gb",
        "date": str(datetime.now()),
    }

    return TraversalResult(
        initfinaldiff=df_summary.to_dict(),
        final_concs=mean_concs.to_dict(),
        data=samples,
        metadata=metadata,
    )

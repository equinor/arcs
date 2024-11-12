from __future__ import annotations
from dataclasses import dataclass
from typing import Literal, Any

from chempy.equilibria import Equilibrium, EqSystem
from chempy import Substance
import copy
import networkx as nx
import itertools as it

from tqdm import tqdm
import platform
import psutil

from datetime import datetime
import numpy as np
import pandas as pd

from arcs.model import get_graph, get_reactions


@dataclass
class TraversalResult:
    initfinaldiff: Any
    final_concs: Any
    metadata: dict[str, Any]
    data: dict[str, Any]


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
    concs = copy.deepcopy(init_concs)  # don't modify the original
    if co2 is False and "CO2" in concs:
        del concs["CO2"]  # CO2 will always be too large as it is the background
    # house keeping:
    num_not_zero = len([x for x in concs.values() if x > 0])
    if max_compounds > num_not_zero:
        max_compounds = num_not_zero

    # scale the probabilities accordingly based upon a ceiling percentage
    median_conc = np.median(
        [v for v in concs.values() if v > 0]
    )  # median > mean for this
    # new_concs = {}
    above_ceiling = {
        k: v for k, v in concs.items() if v > (median_conc * (1 + (ceiling / 100)))
    }
    # modify the ceiling by scaling it down to a suitable value
    # should still max out if concentrations become way to high
    for k, v in above_ceiling.items():
        concs[k] = v * scale_highest

    # get the probabilities based upon relative concentrations:
    p_1 = {k: v / sum(concs.values()) for k, v in concs.items()}
    # now filter based upon the probability threshold:
    p_2 = {k: v for k, v in p_1.items() if v > probability_threshold}
    p_3 = {k: v / sum(p_2.values()) for k, v in p_2.items()}
    # make a list of choices based upon the probabilities
    available = list(
        rng.choice(list(p_3.keys()), 100, p=list(p_3.values()))
    )  # make this list length of the nodes
    # now make a list max_compounds long of random choices based on available
    choices = {}
    for c in range(max_compounds):
        if c == 0:
            c1 = rng.choice(available)
            choices[c1] = p_3[c1]
        else:
            try:
                for i in range(available.count(list(choices)[c - 1])):
                    available.remove(list(choices)[c - 1])
                try:
                    c2 = rng.choice(available)
                    choices[c2] = p_3[c2]
                except Exception:
                    pass
            except Exception:
                pass

    return choices


def _length_multiplier(candidate, *, rank_small_reactions_higher: bool):
    if rank_small_reactions_higher:
        return len(list(candidate))
    else:
        return 1


def _get_weighted_reaction_rankings(
    tempreature: int,
    pressure: int,
    choices: list[str],
    *,
    max_rank: int,
    method: Literal["Bellman-Ford", "Dijkstra"],
    rank_small_reactions_higher: bool,
):
    graph = get_graph(tempreature, pressure)
    rankings = {}
    if len(choices) > 1:
        possibilities = list(
            nx.shortest_paths.all_shortest_paths(
                graph, list(choices)[0], list(choices)[1], method=method
            )
        )

        for x in possibilities:
            candidates = list(graph[x[1]])
            if len(choices) > 2:
                for c in list(choices)[2:]:
                    if c in candidates:
                        weight = graph.get_edge_data(x[0], x[1])[0][
                            "weight"
                        ] * 10 ** _length_multiplier(
                            graph[x[1]],
                            rank_small_reactions_higher=rank_small_reactions_higher,
                        )
                        rankings[x[1]] = {
                            "candidates": candidates,
                            "weight": weight,
                        }
            else:
                weight = graph.get_edge_data(x[0], x[1])[0][
                    "weight"
                ] * 10 ** _length_multiplier(
                    graph[x[1]],
                    rank_small_reactions_higher=rank_small_reactions_higher,
                )
                rankings[x[1]] = {"candidates": candidates, "weight": weight}
    if rankings:
        sorted_rankings = (
            pd.DataFrame(rankings).sort_values(by="weight", axis=1).to_dict()
        )
        topranks = [
            x for i, x in enumerate(sorted_rankings) if i <= max_rank
        ]  # need to sort first
        rankings = {x: rankings[x] for x in topranks}
        return rankings
    else:
        return None


def _generate_eqsystem(index: int, temperature: int, pressure: int) -> EqSystem | None:
    charged_species = {
        "CO3H": -1,
        "NH4": +1,
        "NH2CO2": -1,
    }  # this needs to be added to the arguments
    rs = get_reactions(temperature, pressure)[index]
    r = rs["e"].reac
    p = rs["e"].prod
    k = rs["k"]
    substances = {}
    for n in list(it.chain(*[list(r) + list(p)])):
        if n in list(charged_species.keys()):
            s = Substance.from_formula(n, **{"charge": charged_species[n]})
            substances[s.name] = s
        else:
            s = Substance.from_formula(
                n, **{"charge": 0}
            )  # ,charge=0) # charge buggers up everything have removed for now....
            substances[s.name] = s
    eql = Equilibrium(reac=r, prod=p, param=k)
    try:
        return EqSystem([eql], substances)  # might not just be able to try a return...
    except Exception:
        return None


def _equilibrium_concentrations(
    concs: dict[str, float], eq: EqSystem
) -> tuple[dict[str, float], str]:
    # something is going wrong here...
    fc = copy.deepcopy(concs)
    try:
        x, sol, sane = eq.root(fc)
        assert sol["success"] and sane
        for n, c in enumerate(x):
            fc[eq.substance_names()[n]] = c

        concs = fc
        eq = eq.string()
    except Exception:
        concs = fc
        eq = None
    return (concs, eq)


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
):
    final_concs = {0: copy.deepcopy(concs)}
    reactionstats = {0: None}

    for ip in range(1, path_depth + 1):
        fcs = copy.deepcopy(final_concs[ip - 1])
        try:
            choices = _get_weighted_random_compounds(
                temperature,
                pressure,
                fcs,
                max_compounds=max_compounds,
                probability_threshold=probability_threshold,
                co2=co2,
                scale_highest=scale_highest,
                ceiling=ceiling,
                rng=rng,
            )
        except Exception:
            path_depth = ip + 1
            break
        if len(choices) <= 1:  # not sure this is necessary....
            path_depth = ip + 1
            break
        rankings = _get_weighted_reaction_rankings(
            temperature,
            pressure,
            choices,
            max_rank=max_rank,
            method=method,
            rank_small_reactions_higher=rank_small_reactions_higher,
        )
        if not rankings:
            break
        weights = {k: 1 / rankings[k]["weight"] for k in rankings}
        probabilities = {k: v / sum(weights.values()) for k, v in weights.items()}
        chosen_reaction = rng.choice(
            [
                rng.choice(
                    list(probabilities.keys()),
                    len(probabilities),
                    p=list(probabilities.values()),
                )
            ][0]
        )

        eqsyst = _generate_eqsystem(chosen_reaction, temperature, pressure)
        # if reaction was previous reaction then break
        path_available = [r for r in reactionstats.values() if r is not None]
        if path_available:
            if (
                eqsyst.string() == path_available[-1]
                and eqsyst.string() == path_available[-1]
            ):
                break  # extra break

        final_concs[ip], reactionstats[ip] = _equilibrium_concentrations(fcs, eqsyst)

    return {
        "data": final_concs[list(final_concs)[-1]],
        "equation_statistics": [r for r in reactionstats.values() if r is not None],
        "path_length": len([r for r in reactionstats.values() if r is not None]),
    }


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
    method: Literal["Bellman-Ford"],
    rng: np.random.Generator,
) -> dict[int, Any]:
    init_concs = copy.deepcopy(concs)
    result_dict = {
        0: {"data": init_concs, "equation_statistics": [], "path_length": None}
    }
    with tqdm(
        total=sample_length,
        bar_format="progress: {desc:<10}|{bar:50}|",
        ascii=" >=",
        position=0,
        leave=False,
    ) as pbar:
        for sample in range(sample_length):
            result_dict[sample + 1] = _random_walk(
                temperature,
                pressure,
                init_concs,
                probability_threshold=probability_threshold,
                path_depth=path_depth,
                max_compounds=max_compounds,
                max_rank=max_rank,
                co2=co2,
                scale_highest=scale_highest,
                ceiling=ceiling,
                method=method,
                rank_small_reactions_higher=rank_small_reactions_higher,
                rng=rng,
            )
            pbar.update(1)
    return result_dict


def traverse(
    temperature: int,
    pressure: int,
    concs: dict[str, float] | None = None,
    *,
    co2: bool = False,
    max_compounds: int = 5,
    probability_threshold: float = 0.05,
    max_rank: int = 5,
    sample_length: int = 1000,
    path_depth: int = 20,
    ceiling: int = 2000,
    scale_highest: float = 0.1,
    rank_small_reactions_higher: bool = True,
    method: Literal["Bellman-Ford", "Dijkstra"] = "Bellman-Ford",
    rng: np.random.Generator | None = None,
) -> TraversalResult:
    if rng is None:
        rng = np.random.default_rng()

    concs = {**concs}
    concstring = pd.Series({k: v for k, v in concs.items() if v > 0}) / 1e-6
    if "CO2" in concstring:
        del concstring["CO2"]

    path_lengths = []
    data = _sample(
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
    )

    reformatted = [{x: v for x, v in data[i]["data"].items()} for i in data]
    mean = (
        pd.Series(
            {k: v for k, v in pd.DataFrame(reformatted).mean().items() if v > 0.5e-6}
        )
        / 1e-6
    )

    final_concs = mean.to_dict()
    diff_concs = pd.Series(mean.to_dict()) - pd.Series(
        {k: v / 1e-6 for k, v in concs.items()}
    )
    ift = pd.DataFrame(
        [
            {k: v / 1e-6 for k, v in concs.items() if v > 0},
            mean.to_dict(),
            diff_concs.to_dict(),
        ],
        index=["initial", "final", "change"],
    ).T
    initfinaldiff = ift.dropna(how="all").fillna(0.0).to_dict()
    avgpathlength = np.median(
        [data[i]["path_length"] for i in data if data[i]["path_length"] is not None]
    )

    path_lengths.append(avgpathlength)

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
        initfinaldiff=initfinaldiff,
        final_concs=final_concs,
        data=data,
        metadata=metadata,
    )

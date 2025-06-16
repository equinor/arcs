from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, TypedDict

from chempy.equilibria import Equilibrium, EqSystem
from chempy import Substance
import copy
import itertools as it
import concurrent.futures
import platform
import psutil

from datetime import datetime
import numpy as np
import pandas as pd

from arcs.model import get_reaction_compounds, get_table, get_reactions, Table
import warnings

warnings.filterwarnings("ignore")


if TYPE_CHECKING:
    import numpy.typing as npt


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
) -> dict[str, float]:
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


def _get_weighted_reaction_rankings(
    choices: list[str],
    *,
    max_rank: int,
    table: Table,
    reaction_compounds: dict[int, set[str]],
) -> tuple[npt.NDArray[np.int16], npt.NDArray[np.float64]] | None:
    reactions, weights = table[(str(choices[0]), str(choices[1]))]

    new_weights: list[np.float64] = []
    new_reacts: list[np.int16] = []

    for index, reaction in enumerate(reactions):
        if len(new_weights) > max_rank:
            break

        if len(choices) <= 2 or any(
            c in reaction_compounds[int(reaction)] for c in choices[2:]
        ):
            new_reacts.append(reaction)
            new_weights.append(weights[index])

    weights = np.array(new_weights)
    reactions = np.array(new_reacts)

    sum = np.sum(weights)
    if not sum:
        return None

    weights /= sum
    return (reactions, weights)


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
    iter: int,
    max_compounds: int,
    max_rank: int,
    co2: bool,
    scale_highest: float,
    ceiling: int,
    rng: np.random.Generator,
    reactions: dict[int, Any],
    table: Table,
    reaction_compounds: dict[int, set[str]],
) -> _RandomWalk:
    conc_history = [concs]
    reaction_history: list[str] = []

    for _ in range(iter):
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
        rankings = _get_weighted_reaction_rankings(
            list(choices),
            max_rank=max_rank,
            table=table,
            reaction_compounds=reaction_compounds,
        )
        if rankings is None:
            continue
        chosen_reaction = rng.choice(
            rng.choice(
                rankings[0],
                len(rankings[0]),
                p=rankings[1],
            )
        )

        eqsyst = _generate_eqsystem(
            int(chosen_reaction), temperature, pressure, reactions=reactions
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
    samples: int,
    iter: int,
    ceiling: int,
    scale_highest: float,
    rng: np.random.Generator,
    reactions: dict[int, Any],
    table: Table,
    nproc: int,
    reaction_compounds: dict[int, set[str]],
) -> dict[int, Any]:
    config = {
        "temperature": temperature,
        "pressure": pressure,
        "concs": concs,
        "co2": co2,
        "max_compounds": max_compounds,
        "probability_threshold": probability_threshold,
        "max_rank": max_rank,
        "iter": iter,
        "ceiling": ceiling,
        "scale_highest": scale_highest,
        "rng": rng,
        "reactions": reactions,
        "table": table,
        "reaction_compounds": reaction_compounds,
    }
    if nproc == 0:
        nproc = psutil.cpu_count() or 1  # because cpu_count can return None

    chunk_size = samples // nproc

    result_dict: dict[int, _RandomWalk] = {
        0: {"data": concs, "equation_statistics": [], "path_length": 0}
    }

    if nproc == 1 or samples < nproc:
        print("Single process")
        result_dict.update(_sample_chunk(0, samples, config))
    else:
        print(f"CPU count: {nproc}, chunk size: {chunk_size}")
        rngs = rng.spawn(nproc)
        with concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as executor:
            futures = {
                executor.submit(
                    _sample_chunk, chunk_id, chunk_size, {**config, "rng": rng}
                ): chunk_id
                for chunk_id, rng in enumerate(rngs)
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
    samples: int = 1000,
    iter: int = 5,
    ceiling: int = 500,
    scale_highest: float = 0.1,
    rank_small_reactions_higher: bool = True,
    rng: np.random.Generator | None = None,
    reactions: dict[int, Any] | None = None,
    table: Table | None = None,
    nproc: int = 0,
) -> TraversalResult:
    if rng is None:
        rng = np.random.default_rng()
    if table is None:
        table = get_table(
            temperature,
            pressure,
            rank_small_reactions_higher=rank_small_reactions_higher,
        )
    if reactions is None:
        reactions = get_reactions(temperature, pressure)

    concstring = pd.Series({k: v for k, v in concs.items() if v > 0})
    if "CO2" in concstring:
        del concstring["CO2"]

    path_lengths = []
    results = _sample(
        temperature,
        pressure,
        concs,
        co2=co2,
        max_compounds=max_compounds,
        probability_threshold=probability_threshold,
        max_rank=max_rank,
        samples=samples,
        iter=iter,
        ceiling=ceiling,
        scale_highest=scale_highest,
        rng=rng,
        reactions=reactions,
        table=table,
        nproc=nproc,
        reaction_compounds=get_reaction_compounds(reactions),
    )

    mean_concs = pd.DataFrame([s["data"] for s in results.values()]).mean()
    df_summary = pd.DataFrame({"initial": concs, "final": mean_concs}) * 1e6
    df_summary = df_summary.dropna(how="all").fillna(0.0)
    df_summary["change"] = df_summary["final"] - df_summary["initial"]

    df_summary = df_summary.loc[(df_summary.abs() >= 1e-6).any(axis=1)]

    avg_path_length = np.median(
        [s["path_length"] for s in results.values() if s["path_length"] is not None]
    )
    path_lengths.append(avg_path_length)

    metadata = {
        "arcs_version": "1.4.0",
        "avg_path_length": np.mean(path_lengths),
        "co2": co2,
        "max_compounds": max_compounds,
        "probability_threshold": probability_threshold,
        "max_rank": max_rank,
        "samples": samples,
        "iter": iter,
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
        data=results,
        metadata=metadata,
    )

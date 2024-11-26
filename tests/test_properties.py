from arcs.traversal import (
    _sample,
    _get_weighted_reaction_rankings,
    _get_weighted_random_compounds,
)
from hypothesis import given, settings
from hypothesis.strategies import integers, floats
import numpy as np
from arcs.model import get_graph, get_reactions
import math
from random import random

@given(integers(min_value=1, max_value=1500))
@settings(max_examples=50)
def test_sample_length(sample_length: int):
    concentrations = {
        "CO2": 1.0,
        "H2O": 30.0e-6,
        "O2": 10.0e-6,
        "SO2": 10.0e-6,
        "NO2": 0,
        "H2S": 10.0e-6,
    }
    temperature = 300
    pressure = 10
    sample = _sample(
        temperature=temperature,
        pressure=pressure,
        concs=concentrations,
        co2=True,
        max_compounds=5,
        probability_threshold=0.05,
        max_rank=5,
        sample_length=sample_length,
        path_depth=5,
        ceiling=2000,
        scale_highest=0.1,
        rank_small_reactions_higher=True,
        method="Bellman-Ford",
        rng=np.random.default_rng([0]),
        reactions=get_reactions(temperature, pressure),
        graph=get_graph(temperature, pressure),
    )
    assert len(sample.keys()) == sample_length + 1


@given(
    floats(min_value=1e-30, max_value=1e30),
    floats(min_value=1e-30, max_value=1e30),
    floats(min_value=1e-30, max_value=1e30),
    floats(min_value=1e-30, max_value=1e30),
)
def test_probability_sums_to_one(
    conc1: float, conc2: float, conc3: float, conc4: float
):
    temperature = 300
    pressure = 10
    concentrations = {"SO2": conc1, "NO2": conc2, "H2S": conc3, "H2O": conc4}

    weighted_random_compounds = _get_weighted_random_compounds(
        temperature=temperature,
        pressure=pressure,
        init_concs=concentrations,
        co2=False,
        max_compounds=5,
        probability_threshold=0.05,
        scale_highest=0.1,
        ceiling=2000,
        rng=np.random.default_rng([0]),
    )

    sum_probabilities = 0
    for compound in weighted_random_compounds.keys():
        sum_probabilities += weighted_random_compounds[compound]
    assert math.isclose(sum_probabilities, 1)


@given(integers(min_value=0))
def test_number_of_rankings_doesnt_exceed_max_rank(max_rank: int):
    temperature = 300
    pressure = 10
    graph = get_graph(temperature, pressure)
    concentrations = {"SO2": 10e-6, "NO2": 50e-6, "H2S": 30e-6, "H2O": 20e-6}

    choices = _get_weighted_random_compounds(
        temperature=temperature,
        pressure=pressure,
        init_concs=concentrations,
        co2=False,
        max_compounds=5,
        probability_threshold=0.05,
        scale_highest=0.1,
        ceiling=2000,
        rng=np.random.default_rng([0]),
    )

    number_of_rankings = len(
        _get_weighted_reaction_rankings(
            tempreature=temperature,
            pressure=pressure,
            choices=choices,
            max_rank=max_rank,
            method="Bellman-Ford",
            rank_small_reactions_higher=True,
            graph=graph,
        )
    )

    assert number_of_rankings <= max_rank + 1

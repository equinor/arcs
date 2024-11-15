import numpy as np
import pandas as pd
from arcs.traversal import traverse
from arcs.analysis import AnalyseSampling
from chempy import Equilibrium
import networkx as nx
import pytest
import math


def test_snapshot(snapshot):
    temperature = 300
    pressure = 10

    concs = {
        "CO2": 1.0,
        "H2O": 30.0e-6,
        "O2": 10.0e-6,
        "SO2": 10.0e-6,
        "NO2": 0,
        "H2S": 10.0e-6,
    }

    results = traverse(
        temperature,
        pressure,
        concs,
        sample_length=100,
        path_depth=5,
        ceiling=500,
        scale_highest=0.1,
        max_rank=5,
        max_compounds=5,
        method="Dijkstra",
        rng=np.random.default_rng([0]),
        probability_threshold=0.05,
    )

    analysis = AnalyseSampling(results.data, markdown=True)
    analysis.reaction_statistics()

    df = pd.DataFrame(analysis.stats)

    snapshot.assert_match(df.to_csv(), "snapshot.csv")


@pytest.mark.parametrize(
    "k,expected_left_shift,expected_right_shift",
    [(1e-12, True, False), (1, False, False), (1e12, False, True)],
)
def test_synthetic(k, expected_left_shift, expected_right_shift):
    reactions = {
        1: {
            "e": Equilibrium.from_string("H2 + Br2 = 2 HBr"),
            "k": k,
        }
    }
    graph = nx.MultiDiGraph(directed=True)

    for i, reac in reactions.items():
        f_cost = 1
        b_cost = 1
        reactants = list(reac["e"].reac)
        products = list(reac["e"].prod)
        graph.add_weighted_edges_from(
            [c, i, f_cost] for c in reactants
        )  # reactants -> reaction
        graph.add_weighted_edges_from(
            [i, c, b_cost] for c in reactants
        )  # reaction -> reactants
        graph.add_weighted_edges_from(
            [i, c, f_cost] for c in products
        )  # reaction -> products
        graph.add_weighted_edges_from(
            [c, i, b_cost] for c in products
        )  # products -> reaction

    concs = {
        "H2": 1e-6,
        "Br2": 1e-6,
        "HBr": 1e-6,
    }

    res = traverse(
        temperature=10,
        pressure=10,
        concs=concs,
        graph=graph,
        reactions=reactions,
        sample_length=1,
        path_depth=1,
    )

    first_product_change = res.initfinaldiff["change"][products[0]]
    if not expected_left_shift and not expected_right_shift:
        assert math.isclose(first_product_change, 0, abs_tol=1e-15)

    if expected_left_shift and not expected_right_shift:
        assert first_product_change < 0

    if not expected_left_shift and expected_right_shift:
        assert first_product_change > 0

import numpy as np
import pytest

from arcs.analysis import AnalyseSampling
from arcs.traversal import traverse


@pytest.fixture
def analysis():
    data = {
        0: {
            "data": {"A": 1.0 * 1e-6, "B": 2.0 * 1e-6, "C": 3.0 * 1e-6},
            "equation_statistics": ["A + B = C; 1.0", "B + C = D; 2.0"],
        },
        1: {
            "data": {"A": 1.5 * 1e-6, "B": 2.5 * 1e-6, "C": 3.5 * 1e-6},
            "equation_statistics": ["A + B = C; 1.0", "C + D = E; 3.0"],
        },
        2: {
            "data": {"A": 2.0 * 1e-6, "B": 3.0 * 1e-6, "C": 4.0 * 1e-6},
            "equation_statistics": ["A + B = C; 1.0", "B + C = D; 2.0"],
        },
    }
    return AnalyseSampling(data)


def test_mean_sampling(analysis):
    analysis.mean_sampling()

    expected_final_concs = {"A": 1.5, "B": 2.5, "C": 3.5}
    expected_mean_data = {
        "A": {"value": 0.5, "variance": 0.25},
        "B": {"value": 0.5, "variance": 0.25},
        "C": {"value": 0.5, "variance": 0.25},
    }

    assert analysis.final_concs == pytest.approx(expected_final_concs, rel=1e-6)
    assert analysis.mean_data["A"] == pytest.approx(expected_mean_data["A"], rel=1e-6)
    assert analysis.mean_data["B"] == pytest.approx(expected_mean_data["B"], rel=1e-6)
    assert analysis.mean_data["C"] == pytest.approx(expected_mean_data["C"], rel=1e-6)


def test_get_stats(analysis):
    equations = [
        ["A + B = C; k=1.0", "B + C = D; k=2.0"],
        ["A + B = C; k=1.0", "C + D = E; k=3.0"],
        ["A + B = C; k=1.0", "B + C = D; k=2.0"],
    ]

    expected_stats = {
        "index": {0: "A + B = C", 1: "B + C = D", 2: "C + D = E"},
        "k": {0: " k=1.0", 1: " k=2.0", 2: " k=3.0"},
        "frequency": {0: 3, 1: 2, 2: 1},
    }

    stats = analysis._get_stats(equations)
    assert stats == expected_stats


def test_reaction_paths(analysis):
    analysis.reaction_paths()

    expected_common_paths = {
        "paths": {0: "A + B = C  \n B + C = D "},
        "k": {0: "1.0 \n 2.0"},
        "frequency": {0: 1},
    }

    assert analysis.common_paths == expected_common_paths


def test_concentration_does_not_change():
    # Given the configuration and concentration there should not be any reactions hence no change in concentratons

    temperature = 280
    pressure = 26

    concs = {
        "CO2": 0.0,
        "H2O": 6e-6,
        "O2": 0,
        "SO2": 10.0e-6,
        "NO2": 0,
        "H2S": 10.0e-6,
    }

    results = traverse(
        temperature,
        pressure,
        concs,
        samples=100,
        iter=5,
        ceiling=500,
        scale_highest=0.1,
        max_rank=5,
        max_compounds=5,
        rng=np.random.default_rng([0]),
        probability_threshold=0.05,
        nproc=1,
    )

    conc_in = results.initfinaldiff["initial"]
    conc_out = results.initfinaldiff["final"]

    assert all(np.isclose(list(conc_in.values()), list(conc_out.values())))

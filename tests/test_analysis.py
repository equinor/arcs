import pytest

from arcs.analysis import AnalyseSampling


@pytest.fixture
def analysis():
    data = {
        0: {
            "data": {"A": 1.0, "B": 2.0, "C": 3.0},
            "equation_statistics": ["A + B = C; 1.0", "B + C = D; 2.0"],
        },
        1: {
            "data": {"A": 1.5, "B": 2.5, "C": 3.5},
            "equation_statistics": ["A + B = C; 1.0", "C + D = E; 3.0"],
        },
        2: {
            "data": {"A": 2.0, "B": 3.0, "C": 4.0},
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

    assert analysis.final_concs == expected_final_concs
    assert analysis.mean_data == expected_mean_data


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

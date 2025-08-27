from arcs.traversal import traverse
import numpy as np


def test_balance_nitrogen():
    """This is the specific case were the bug was noticed"""

    concs = {
        "O2": 10.0e-6,
        "H2O": 30.0e-6,
        "SO2": 10.0e-6,
        "NO2": 2.0e-6,
        "H2S": 1.0e-6,
    }

    results = traverse(
        temperature=300,
        pressure=10,
        concs=concs,
        samples=10,
        iter=5,
        ceiling=500,
        scale_highest=0.1,
        max_rank=5,
        max_compounds=5,
        rng=np.random.default_rng([0]),
        probability_threshold=0.05,
        nproc=1,
    ).final_concs

    substance_multipliers = {
        "HNO2": 1,
        "HNO3": 1,
        "N2": 2,
        "NO2": 1,
        "NOHSO4": 1,
    }

    actual = sum(
        substance_multipliers.get(subst, 0) * amount
        for subst, amount in results.items()
    )

    assert np.isclose(actual, 2.0e-6)

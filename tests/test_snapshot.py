import numpy as np
import pandas as pd
from arcs.traversal import traverse
from arcs.analysis import AnalyseSampling


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

import pandas as pd
from arcs.traversal import traverse
from arcs.analysis import AnalyseSampling
import warnings
import argparse
import time

warnings.filterwarnings("ignore")


def run_simulation(samples, iter, nproc):
    start = time.time()
    temperature = 300
    pressure = 10
    concs = {"SO2": 10e-6, "NO2": 0, "H2S": 10e-6, "H2O": 30e-6, "O2": 10e-6}

    results = traverse(
        temperature,
        pressure,
        concs,
        samples=samples,
        nproc=nproc,
        iter=iter,
    )

    analysis = AnalyseSampling(results.data)
    analysis.reaction_statistics()
    df = pd.DataFrame(analysis.stats)
    print(df)

    df_d = pd.DataFrame(results.initfinaldiff)
    print(df_d)

    end = time.time()
    print(f"Simulation took {end - start:.2f} seconds.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--iter", type=int, default=5)
    parser.add_argument("--nproc", type=int, default=0)

    args = parser.parse_args()
    run_simulation(args.samples, args.iter, args.nproc)

import pandas as pd
from arcs.traversal import traverse
import warnings
import argparse
import time

warnings.filterwarnings("ignore")


def run_simulation(sample_length=100):
    start = time.time()
    temperature = 300
    pressure = 10
    concs = {"SO2": 10e-6, "NO2": 0, "H2S": 10e-6, "H2O": 30e-6, "O2": 10e-6}

    results = traverse(
        temperature,
        pressure,
        concs,
        sample_length=sample_length,
        nproc=0,
    )

    df_d = pd.DataFrame(results.initfinaldiff)
    print(df_d)

    end = time.time()
    print(f"Simulation took {end - start:.2f} seconds.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample_length", type=int, default=100)
    args = parser.parse_args()
    run_simulation(args.sample_length)

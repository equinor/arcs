import numpy as np
from numpy.random import PCG64
from arcs.traversal import (
    traverse,
)
import pytest
import time


@pytest.mark.parametrize(
    "sample_length, number_of_runs",
    [(10, 100), (20, 100), (40, 100), (80, 100), (160, 100)],
)
def test_sample_length(sample_length, number_of_runs):
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

    final_concentrations =[]
    start_time = time.time()
    for i in range(number_of_runs):
        results = traverse(
            temperature,
            pressure,
            concs,
            sample_length=sample_length,
            path_depth=5,
            ceiling=500,
            scale_highest=0.1,
            max_rank=5,
            max_compounds=5,
            method="Dijkstra",
            rng=np.random.Generator(PCG64()),
            probability_threshold=0.05,
        )
        #print(f"Final concentrations are: {results.final_concs}")
        final_concentrations.append(results.final_concs)

    run_time = time.time() - start_time
    compound_values = {}

    # Populate the compound_values with results from each run
    for run in final_concentrations:
        for compound, value in run.items():
            if compound not in compound_values:
                compound_values[compound] = []
            compound_values[compound].append(value)

    print(f"\n This is the run with  the sample length of {sample_length} and it was run {number_of_runs} times.\n")
    print(f"The total run time is {run_time}.\n")
    # Calculate mean and standard deviation for each compound
    for compound, values in compound_values.items():
        mean_value = np.mean(values)
        std_deviation = np.std(values)
        print(f"{compound} - Mean: {mean_value}, Std: {std_deviation}")

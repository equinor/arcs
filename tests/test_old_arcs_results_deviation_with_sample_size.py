import pickle

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from arcs.traversal import Traversal

from arcs.setup_functions import GenerateInitialConcentrations

from arcs.analysis import AnalyseSampling

sample_sizes = [10, 100, 1000]

def test_mean_std_across_arcs_sampling():
    pressure = 10
    temperature = 300
    concentrations = {
        "CO2": 1.0,
        "H2O": 500.0e-6,
        "SO2": 100.0e-6,
        "NO2": 30.0e-6,
        "H2S": 50.0e-6,
        "O2": 100.0e-6,
        "H2SO4": 100.0e-6,
    }

    mean_values = {}
    std_values = {}
    file_location = "../app/data/"

    for s in sample_sizes:
        g = pickle.load(open(file_location + "SCAN_graph.p", "rb"))
        t = Traversal(graph=g, reactions=file_location + "SCAN_reactions.p")

        gic = GenerateInitialConcentrations(g)
        gic.all_zero(include_co2=False)
        gic.update_ic(concentrations)
        concs = gic.ic
        settings = {
            "nprocs": 1,
            "sample_length": s,
            "max_rank": 10,
            "max_compounds": 5,
            "probability_threshold": 0.1,
            "path_depth": 5,
            "ceiling": 2000,
            "scale_highest": 0.2,
            "rank_small_reactions_higher": True
        }
        results = t.run([300], [10], ic=concs, **settings, )

        # results = traverse(
        #     temperature,
        #     pressure,
        #     concentrations,
        #     samples=s,
        #     iter=100,
        # )
        analysis = AnalyseSampling(t.data)
        analysis.output_concentrations_std_mean()
        std = analysis.final_concs_std
        mean = analysis.final_concs_mean
        print(s, std, mean)

        for component in mean.keys():
            if component not in mean_values:
                mean_values[component] = []
                std_values[component] = []

            mean_values[component].append(mean[component])
            std_values[component].append(std[component])

    component_list = ["H2O", "H2SO4"]
    plot_mean_std_from_component_list(component_list, mean_values, std_values)


def plot_mean_std_from_component_list(component_list, mean_values, std_values):
    # Plotting
    plt.figure(figsize=(12, 6))

    # Check sizes before plotting
    for component in component_list:
        print(f"Component: {component}")
        print(f"Sample Sizes: {sample_sizes}")
        print(f"Mean Values: {mean_values[component]}")
        print(f"Std Values: {std_values[component]}")

        # Check if lengths match
        if len(sample_sizes) == len(mean_values[component]) and len(sample_sizes) == len(std_values[component]):
            plt.errorbar(sample_sizes, mean_values[component], yerr=std_values[component],
                         label=component, capsize=5, fmt='o', markersize=5)
        else:
            print(f"Warning: Sizes do not match for component {component}. Skipping this component.")

    # Adding labels and title
    plt.xscale('log')
    plt.xlabel('Sample Size (s)')
    plt.ylabel('Concentration')
    plt.title('Mean and Standard Deviation of Concentrations vs Sample Size')
    plt.legend()
    plt.grid()
    plt.show()







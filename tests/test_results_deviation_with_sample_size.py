import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from arcs.traversal import traverse
from arcs.analysis import AnalyseSampling

sample_sizes = [100, 250, 500, 1000, 2000]

def test_mean_std_across_arcs_sampling():
    pressure = 10
    temperature = 300
    concentrations = {
        "CO2": 1.0,
        "H2O": 40.0e-6,
        "O2": 10.0e-6,
        "SO2": 10.0e-6,
        "NO2": 5.0e-6,
        "H2S": 5.0e-6,
    }

    mean_values = {}
    std_values = {}

    for s in sample_sizes:
        results = traverse(
            temperature,
            pressure,
            concentrations,
            samples=s,
        )
        analysis = AnalyseSampling(results.data)
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


def test_std_of_final_concentrations_when_running_arcs_several_times(repeats=10):
    pressure = 10
    temperature = 300
    concentrations = {
        "CO2": 1.0,
        "H2O": 40.0e-6,
        "O2": 10.0e-6,
        "SO2": 10.0e-6,
        "NO2": 5.0e-6,
        "H2S": 5.0e-6,
    }
    mean_values = {}
    std_values = {}

    for s in sample_sizes:
        s = int(s/repeats)
        final_conc_n_runs = {}
        for n in range(repeats):
            results = traverse(
                temperature,
                pressure,
                concentrations,
                samples=s,
            )
            for chemical in results.final_concs.keys():
                if chemical not in final_conc_n_runs:
                    final_conc_n_runs[chemical] = []
                final_conc_n_runs[chemical].append(results.final_concs[chemical]*1e6)

        for chemical in final_conc_n_runs.keys():
            if chemical not in mean_values:
                mean_values[chemical] = []
                std_values[chemical] = []
            mean_values[chemical].append(np.mean(final_conc_n_runs[chemical]))
            std_values[chemical].append(np.std(final_conc_n_runs[chemical]))

    plot_mean_std_from_component_list(["H2O", "H2SO4"], mean_values, std_values)











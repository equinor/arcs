import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from arcs.traversal import traverse
from arcs.analysis import AnalyseSampling

sample_sizes = [10, 100, 500]

def test_sample_size_deviation():
    #sample_sizes = [500, 1000, 1500, 2000]
    #n = 100

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
    plt.xlabel('Sample Size (s)')
    plt.ylabel('Concentration')
    plt.title('Mean and Standard Deviation of Concentrations vs Sample Size')
    plt.legend()
    plt.grid()
    plt.show()


def test_std_of_final_concentrations_when_running_arcs_several_times(repeats=100):
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
                final_conc_n_runs[chemical].append(results.final_concs[chemical])

        for chemical in final_conc_n_runs.keys():
            if chemical not in mean_values:
                mean_values[chemical] = []
                std_values[chemical] = []
            mean_values[chemical].append(np.mean(final_conc_n_runs[chemical]))
            std_values[chemical].append(np.std(final_conc_n_runs[chemical]))

    plot_mean_std_from_component_list(["H2O", "H2SO4"], mean_values, std_values)


def calculate_standard_deviation(data_dict):
    output_concentrations = {}
    for item in data_dict:
        data = data_dict[item]["data"]
        for compound, concentration in data.items():
            if compound not in output_concentrations:
                output_concentrations[compound] = []
            output_concentrations[compound].append(concentration)

    conc_std = {}
    for compound, concentrations in output_concentrations.items():
        conc_std[compound] = float(np.std(concentrations))

    return conc_std


def test_plot_mean_std(result_stats, compounds):
    # Create the plot
    plt.figure(figsize=(8, 6))
    # Collect mean and std for the specified compound
    for c in compounds:
        means = []
        stds = []
        for size in sample_sizes:
            if c in result_stats[size]:
                mean, std = result_stats[size][c]
                means.append(mean*1e6)
                stds.append(std*1e6)
            else:
                raise ValueError(f"{c} not found in sample size {size}")

        plt.errorbar(sample_sizes, means, yerr=stds, label=c, fmt='o-', capsize=5)

    # Add labels and title
    plt.xlabel('Sample Size')
    plt.ylabel('Mean Concentration')
    plt.title(f'Mean Concentration of {compounds} vs Sample Size with Standard Deviation')
    plt.legend()
    plt.grid(True)
    plt.xscale('linear')
    #plt.yscale('log')  # Use logarithmic scale if values vary greatly
    plt.show()


def test_deviations_inside_arcs_with_multiple_runs(n=10):  # Added n as a parameter with a default value
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
        all_means = []
        all_stds = []

        for _ in range(n):  # Run traverse n times for each sample size
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

            all_means.append(mean)
            all_stds.append(std)

        #######
        # Collect all unique components from means and stds
        all_components = set()
        for mean in all_means:
            all_components.update(mean.keys())
        for std in all_stds:
            all_components.update(std.keys())

        # Average the means and standard deviations over n runs for each unique component
        avg_mean = {component: sum(d.get(component, 0) for d in all_means) / n for component in all_components}
        avg_std = {component: sum(d.get(component, 0) for d in all_stds) / n for component in all_components}

        #######

        print(s, avg_std, avg_mean)

        for component in avg_mean.keys():
            if component not in mean_values:
                mean_values[component] = []
                std_values[component] = []

            mean_values[component].append(avg_mean[component])
            std_values[component].append(avg_std[component])

    component_list = ["H2O", "H2SO4"]
    plot_mean_std_from_component_list(component_list, mean_values, std_values)








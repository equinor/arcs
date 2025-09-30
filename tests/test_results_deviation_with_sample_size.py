import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from arcs.traversal import traverse
from arcs.analysis import AnalyseSampling

sample_sizes = [1000]

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

    for s in sample_sizes:
        results = traverse(
            temperature,
            pressure,
            concentrations,
            samples=s,
            iter=100,
        )
        std = results.std_concs
        mean = results.final_concs
        print(s, std, mean)

        for component in mean.keys():
            if component not in mean_values:
                mean_values[component] = []
                std_values[component] = []

            mean_values[component].append(mean[component]*1e6)
            std_values[component].append(std[component]*1e6)

    component_list = ["H2O", "H2SO4", "HNO3"]
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


def merge_two_data_dicts(dict1, dict2):
    # List of dictionaries to merge
    dicts = [dict1, dict2]

    # Merging dictionaries
    merged_dict = {}
    current_index = 0

    for d in dicts:
        for key, value in d.items():
            merged_dict[current_index] = value
            current_index += 1

    return merged_dict

def plot_histogram(df, chemical_name):
    if chemical_name in df.columns:
        plt.figure(figsize=(10, 6))
        plt.hist(df[chemical_name], bins=10, color='blue', alpha=0.7)
        # Calculate the mean
        mean_value = df[chemical_name].mean()

        # Add a line for the mean
        plt.axvline(mean_value, color='red', linestyle='dashed', linewidth=2, label=f'Mean: {mean_value:.2f}')

        plt.title(f'Histogram of {chemical_name}')
        plt.xlabel(chemical_name)
        plt.ylabel('Frequency')
        plt.grid(axis='y', alpha=0.75)
        plt.show()
    else:
        print(f"Chemical '{chemical_name}' not found in the DataFrame.")


def test_distributions_of_final_concentrations():
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
    for s in sample_sizes:
        s = 10
        results = traverse(
            temperature,
            pressure,
            concentrations,
            samples=s,
        )
        data = (pd.DataFrame([s["data"] for s in results.data.values() if s["path_length"] > 0]) * 1e6).fillna(0.0)
        means = data.mean()
        plot_histogram(data, 'H2O')
        plot_histogram(data, 'H2SO4')
        plot_histogram(data, 'HNO3')


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
        for n in range(repeats):
            results = traverse(
                temperature,
                pressure,
                concentrations,
                samples=s,
            )
            chemicals_concentrations = {}
            for chemical in results.final_concs.keys():
                if chemical not in chemicals_concentrations:
                    chemicals_concentrations[chemical] = []
                chemicals_concentrations[chemical].append(results.final_concs[chemical]*1e6)

        for chemical in chemicals_concentrations.keys():
            if chemical not in mean_values:
                mean_values[chemical] = []
                std_values[chemical] = []
            mean_values[chemical].append(np.mean(chemicals_concentrations[chemical]))
            std_values[chemical].append(np.std(chemicals_concentrations[chemical]))

    # plot_mean_std_from_component_list(["H2O", "H2SO4", "HNO3"], mean_values, std_values)
    plot_chemical_data("H2SO4", mean_values)
    plot_chemical_data("H2O", mean_values)


def plot_chemical_data(chemical_name, dataframes):
    """
    Plots the data for the chosen chemical across multiple DataFrames.

    Parameters:
    - chemical_name: str, the name of the chemical to plot.
    - dataframes: list of pd.DataFrame, each DataFrame containing data for plotting.

    Returns:
    - None
    """
    # Create a color map for different datasets
    colors = plt.cm.viridis(np.linspace(0, 1, len(dataframes)))

    # Create a new figure
    plt.figure(figsize=(10, 6))

    # Loop through each DataFrame and plot the data for the chosen chemical
    for i, df in enumerate(dataframes):
        if chemical_name in df.columns:
            plt.hist(df[chemical_name], label=f'Dataset with  {len(df[chemical_name])} samples', color=colors[i])
        else:
            print(f"Warning: {chemical_name} not found in Dataset {i + 1}")

    # Add labels and title
    plt.xlabel('Concentration')
    plt.ylabel('Count')
    plt.title(f'Plot of {chemical_name}')
    plt.legend()
    plt.grid()
    plt.show()









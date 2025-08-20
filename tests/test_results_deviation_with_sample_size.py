import numpy as np

from arcs.traversal import traverse

def test_sample_size_deviation():
    #sample_sizes = [500, 1000, 1500, 2000]
    #n = 100
    sample_sizes = [10, 100, 500]

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

    result_stats = {}

    for s in sample_sizes:
        result_stats[s] = run_arcs(pressure, temperature, concentrations, s)

    print(result_stats)
    test_plot_mean_std(result_stats, ["H2O", "H2SO4"])


def run_arcs(pressure, temperature, concentrations, sample_size):
    results = traverse(
        temperature,
        pressure,
        concentrations,
        samples=sample_size,
    )
    std_mean_per_chem = {}
    std_dict = calculate_standard_deviation(results.data)
    for chem in std_dict.keys():
        std_mean_per_chem[chem] = [results.final_concs[chem], std_dict[chem]]

    return std_mean_per_chem


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


import matplotlib.pyplot as plt

def test_plot_mean_std(result_stats, compounds):
    sample_sizes = list(result_stats.keys())
    means = []
    stds = []
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








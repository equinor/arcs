import numpy as np
import pandas as pd
import json
from arcs.traversal import (
    traverse,
    _get_weighted_reaction_rankings,
    _get_weighted_random_compounds,
    _random_walk,
)
from arcs.model import Table, get_reaction_compounds, get_reactions, get_table
from arcs.analysis import AnalyseSampling
from chempy import Equilibrium
import pytest
import math

from scipy.interpolate import LinearNDInterpolator

pressure_list = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100, 125, 150, 175, 200, 225, 250, 275, 300]
temperature_list = [200, 250, 300, 350, 400]


def test_all_reactions_are_the_same():
    prev_reactions = get_reactions(200, 1)
    for p in pressure_list:
        for t in temperature_list:
            next_reactions = get_reactions(t, p)
            for id in range(len(prev_reactions)):
                assert prev_reactions[id]["e"].prod == next_reactions[id]["e"].prod, (
                    f"Discrepancy found!\n{p=} {t=}\n"
                    f"{prev_reactions[id]["e"].prod=} {next_reactions[id]["e"].prod=}"
                )
                assert prev_reactions[id]["e"].reac == next_reactions[id]["e"].reac, (
                    f"Discrepancy found!\n{p=} {t=}\n"
                    f"{prev_reactions[id]["e"].reac=} {next_reactions[id]["e"].reac=}"
                )
            print(f"All good for {p=} {t=}")
            prev_reactions = next_reactions


def test_gibbs():
    results = []

    for p in pressure_list:
        for t in temperature_list:
            reactions = get_reactions(t, p)
            for reaction_id in range(len(reactions)):
                g_value = reactions[reaction_id]['g']
                results.append({
                    'Reaction_id': reaction_id,
                    'Pressure': p,
                    'Temperature': t,
                    'g': g_value
                })

    # Convert results to a DataFrame for better readability
    df_results = pd.DataFrame(results)

    # Print the DataFrame
    print(df_results)


def test_pressure_and_temperature_values_are_between():
    pressure = 2
    temperature = 250
    pressure_boundaries = [1, 5]
    temperature_boundaries = [200, 300]
    true_reactions = get_reactions(temperature, pressure)

    reactions_0 = get_reactions(200, 1)
    reactions_1 = get_reactions(200, 5)
    reactions_2 = get_reactions(300, 1)
    reactions_3 = get_reactions(300, 5)

    for reaction_id in range(len(true_reactions)):
        assert reactions_0[reaction_id]['g'] <= true_reactions[reaction_id]['g'] <= reactions_2[reaction_id]['g']
        assert reactions_1[reaction_id]['g'] <= true_reactions[reaction_id]['g'] <= reactions_3[reaction_id]['g']


def test_interpolation_gibbs():
    results = []
    pressure_list = [30, 35]
    temperature_list = [200, 250]
    tp_combinations = []

    for p in pressure_list:
        for t in temperature_list:
            tp_combinations.append([t, p])
            reactions = get_reactions(t, p)
            for reaction_id in range(len(reactions)):
                 g_value = reactions[reaction_id]['g']
                 results.append({
                     f"Reaction_id": reaction_id,
                     f"Gibbs_{t}_{p}": g_value
                 })

    # Convert results to a DataFrame for better readability
    df_results = pd.DataFrame(results)

    # Group by 'Reaction_id' and aggregate the Gibbs values
    df_grouped = df_results.groupby("Reaction_id").agg(lambda x: x.dropna().iloc[0] if not x.dropna().empty else None).reset_index()

    # Create a list to hold interpolators for each set of values
    interpolators = [LinearNDInterpolator(tp_combinations, np.array(df_grouped.iloc[i, 1:5])) for i in range(len(df_grouped))]

    # Specify a point (x, y) where you want to get the interpolated values
    target_pressure = 31
    target_temperature = 225

    # Get the interpolated values for the specified point for each set of values
    interpolated_values = [interpolator(target_temperature, target_pressure) for interpolator in interpolators]

    # Print the results
    for j, value in enumerate(interpolated_values):
        print(f"Interpolated value at point ({target_temperature}, {target_pressure}) for Set {j + 1}: {value:.2f}")

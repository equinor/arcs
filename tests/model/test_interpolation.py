import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator

from arcs.model import get_reactions

ALLOWED_PRESSURES_LIST = [1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100, 125, 150, 175, 200, 225, 250, 275, 300]
ALLOWED_TEMPERATURES_LIST = [200, 250, 300, 350, 400]


def test_all_reactions_are_the_same():
    prev_reactions = get_reactions(200, 1)
    for p in ALLOWED_PRESSURES_LIST:
        for t in ALLOWED_TEMPERATURES_LIST:
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

    for p in ALLOWED_PRESSURES_LIST:
        for t in ALLOWED_TEMPERATURES_LIST:
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
    pressure_list = [30, 40]
    temperature_list = [200, 300]
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
    target_pressure = 35
    target_temperature = 250

    # Get the interpolated values for the specified point for each set of values
    linearND_values = [interpolator(target_temperature, target_pressure) for interpolator in interpolators]

    low_pressure_gibbs = [np.interp(target_temperature, temperature_list, [float(df_grouped.iloc[i, 1]), float(df_grouped.iloc[i, 2])]) for i in range(len(df_grouped))]
    high_pressure_gibbs = [np.interp(target_temperature, temperature_list, [float(df_grouped.iloc[i, 3]), float(df_grouped.iloc[i, 4])]) for i in range(len(df_grouped))]

    final_gibbs_linear_interp = [np.interp(target_pressure, pressure_list, [float(low_pressure_gibbs[i]), float(high_pressure_gibbs[i])]) for i in range(len(df_grouped))]

    final_gibbs_log = [log_pressure_interpolation(target_pressure, pressure_list, [float(low_pressure_gibbs[i]), float(high_pressure_gibbs[i])]) for i in range(len(df_grouped))]

    if (np.isclose(linearND_values[i], final_gibbs_linear_interp[i]) for i in range(len(linearND_values))):
        pass
    else:
        print(f"linearND_values = {linearND_values}, final_gibbs = {final_gibbs_linear_interp}")

    if (np.isclose(final_gibbs_linear_interp[i], final_gibbs_log[i]) for i in range(len(linearND_values))):
        pass
    else:
        print(f"final_gibbs_linear_interp = {final_gibbs_linear_interp}, final_gibbs_log = {final_gibbs_log}")

    # Print the results
    #for j, value in enumerate(linearND_values):
    #    print(f"Linearly interpolated value at point ({target_temperature}, {target_pressure}) for Set {j + 1}: {value:.2f}")

    true_reactions = get_reactions(target_temperature, target_pressure)
    for reaction_id in range(len(true_reactions)):
        g_value = true_reactions[reaction_id]['g']
        assert np.isclose(g_value, linearND_values[reaction_id], atol=0.1, rtol=0.01)
        assert np.isclose(g_value, final_gibbs_linear_interp[reaction_id], atol=0.1, rtol=0.01)
        assert np.isclose(g_value, final_gibbs_log[reaction_id], atol=0.1, rtol=0.01)



def log_pressure_interpolation(target: float, x: list, y: list):
    '''
        Interpolating y = a + b log x
        Assuming x = [x0, x1] and y = [y0, y1]
    '''
    b = (y[1] - y[0])/(np.log(x[1]) - np.log(x[0]))
    a = y[0] - b * np.log(x[0])
    log_interpolated_value = a + b * np.log(target)
    return log_interpolated_value



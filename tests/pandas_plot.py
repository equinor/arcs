import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data using pandas
# Replace 'simulation_data.csv' with the path to your actual CSV file
csv_file_path = 'simulation_data.csv'
data = pd.read_csv(csv_file_path)

# Plot the means with the standard deviations as error bars
plt.errorbar(data['sample_length'], data['H2SO4_mean'], yerr=data['H2SO4_std'], fmt='o',
             ecolor='red', capthick=2, capsize=5, label='H2SO4 Mean ± Std')
# Plot the means with the standard deviations as error bars
plt.errorbar(data['sample_length']+1, data['H2SO4_mean'], yerr=4.4*np.power(data['sample_length'],-1/2), fmt='o',
             ecolor='red', capthick=2, capsize=5, label='H2SO4 Mean ± Std')


# Adding labels and title
plt.xlabel('Total Run Time (seconds)')
plt.ylabel('H2SO4 Mean')
plt.title('H2SO4 Mean vs Total Run Time with Error Bars for Std')
plt.legend()
# Set x-axis to log scale
plt.xscale('log')

plt.grid(True)

# Show the plot
plt.show()
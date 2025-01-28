import re


def parse_simulation_results(text_file_path, csv_file_path):
    # Compile regex patterns for different pieces of information
    sample_length_pattern = re.compile(r"sample length of (\d+)")
    total_run_time_pattern = re.compile(r"total run time is ([\d.]+)")
    mean_std_pattern = re.compile(r"(\w+) - Mean: ([\d.]+), Std: ([\d.]+)")

    # Initialize a list to hold all the processed simulation data rows
    simulation_data = []

    # Open and read from the text file
    with open(text_file_path, 'r') as file:
        content = file.read()

        # Split the content by runs
        runs = content.strip().split('\n\n\n')

        # Extract data from each run and form a dictionary for each
        for run in runs:
            # Get sample length and total run time
            sample_length = sample_length_pattern.search(run).group(1)
            total_run_time = total_run_time_pattern.search(run).group(1).strip(".")

            # Find all mean and std values
            mean_std_values = mean_std_pattern.findall(run)

            # Form a dictionary for the current run
            run_data = {
                'sample_length': sample_length,
                'total_run_time': total_run_time,
            }

            # Extend `run_data` with the mean/std values for the gases
            for gas, mean, std in mean_std_values:
                run_data[f'{gas}_mean'] = mean
                run_data[f'{gas}_std'] = std

            # Append the dictionary to the simulation data list
            simulation_data.append(run_data)

    # Define the csv header
    header = list(simulation_data[0].keys())

    # Write the CSV file
    with open(csv_file_path, 'w') as csv_file:
        # Write the header
        csv_file.write(','.join(header) + '\n')

        # Write data rows
        for data in simulation_data:
            csv_file.write(','.join(data.values()) + '\n')


# File paths
input_file = 'simulation_results.txt'  # Replace with your file path
output_file = 'simulation_data.csv'  # Replace with your desired output file path

# Call the function to parse and create the CSV
parse_simulation_results(input_file, output_file)
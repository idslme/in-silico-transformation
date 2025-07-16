import os
import sys
import logging
import pandas as pd
from multiprocessing import Pool
import subprocess

# Set up directories and read the input CSV
os.makedirs('/home/mani/pubchem_rxns_2025_processing/rxn_info_batchpro', exist_ok=True)
unprocess_rxns = pd.read_csv('/home/mani/pubchem_rxns_2025_processing/rxn_diff_euq_complete_unprocessed.csv')
unprocess_rxns.rename(columns={'reaction_cid': 'rxn_id', 'reaction_smiles': 'rxn_smiles'}, inplace=True)
# Define the output directory and subset size
output_dir = '/home/mani/pubchem_rxns_2025_processing/rxn_info_batchpro'
subset_size = 500

# Split the dataframe into smaller subsets and save them in different folders
for i in range(0, len(unprocess_rxns), subset_size):
    subset = unprocess_rxns.iloc[i:i + subset_size]
    subset_dir = os.path.join(output_dir, f"subset_{i // subset_size + 1}")
    os.makedirs(subset_dir, exist_ok=True)  # Create directory for each subset
    output_file = os.path.join(subset_dir, "data.csv")
    subset.to_csv(output_file, index=False)
    print(f"Saved {output_file}")

# Now, collect paths of all the generated subset files
target_dir = '/home/mani/pubchem_rxns_2025_processing/rxn_info_batchpro'
sub_dir_paths = []
for sub_dir in os.listdir(target_dir):
    sub_dir_path = os.path.join(target_dir, sub_dir)
    if os.path.isdir(sub_dir_path):  # Ensure we only list directories
        sub_dir_paths.append(sub_dir_path)
       

file_paths = []
for sub_dir_path in sub_dir_paths:
    # Check if 'reaction_info_with_ids.csv' exists in the sub_dir_path
    if 'reaction_info_with_ids.csv' in os.listdir(sub_dir_path):
        print(f"Skipping folder: {sub_dir_path} as it contains reaction_info_with_ids.csv")
        continue  # Skip this folder and move on to the next one
    
    file_names = os.listdir(sub_dir_path)
    for file_name in file_names:
        if file_name.endswith('.csv'):  # Only process .csv files
            file_path = os.path.join(sub_dir_path, file_name)
            file_paths.append(file_path)

# Path to the reaction info generation script
rxn_info_generation = '/home/mani/ManiProject2025/ManiProject2024/pythonscripts-2025/Metabolite_Predicition_ExposureRX/rxn_info_cli.py'

# Iterate over all subset files and process them
for file_path in file_paths:
    input_file = file_path
    max_num_cpus = '16'
    output_dir = os.path.dirname(file_path)
    
    num_reactions = len(pd.read_csv(input_file))
    
    try:
        print(f"Processing file: {input_file}")
        print(f"Number of reactions to be processed: {num_reactions}")
        # Updated to positional arguments
        result = subprocess.run([
            'python', rxn_info_generation,
            input_file,
            max_num_cpus,
            output_dir
        ], capture_output=True, text=True)

        if result.returncode == 0:
            print(f"Successfully processed {input_file}")
            print(f"Output: {result.stdout}")
        else:
            print(f"Error processing {input_file}")
            print(f"Error Output: {result.stderr}")
    except Exception as e:
        print(f"Failed to run the process for {input_file}. Error: {str(e)}")


directory = '/home/mani/pubchem_rxns_2025_processing/rxn_info_batchpro'
sub_directories = os.listdir(directory)

full_paths = []

for sub_dir in sub_directories:
    sub_dir_path = os.path.join(directory, sub_dir)
    full_paths.append(sub_dir_path)


combined_rxn_info = []

for full_path in full_paths:
    if 'reaction_info_with_ids.csv' in os.listdir(full_path):
        file_path = os.path.join(full_path, 'reaction_info_with_ids.csv')
        df = pd.read_csv(file_path)
        combined_rxn_info.append(df)
rxn_info_combined = pd.concat(combined_rxn_info, ignore_index=True)
# Save the combined DataFrame to a CSV file
rxn_info_combined.to_csv('/home/mani/pubchem_rxns_2025_processing/rxn_info_combined.csv', index=False)
print("Combined reaction info saved to /home/mani/pubchem_rxns_2025_processing/rxn_info_combined.csv")

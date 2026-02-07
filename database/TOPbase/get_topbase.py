# Translate the given bash script into Python code
import requests
import os

# Set up directory
dir_path = "./"

# Define the atomic number range
atomic_numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 26]

# Function to create filenames with leading zeros
def create_filename(prefix, nz, ne):
    return f"{dir_path}/{prefix}.{nz:02d}.{ne:02d}.dat"

# To print the first 10 lines of the response text when an error is detected

# Assuming 'response' is the variable holding the HTTP response from the requests.get call
def print_error(response_text):
    print("*** ERROR ***")
    # Split the response text by new line characters to get a list of lines
    error_lines = response_text.split('\n')
    # Print the first 10 lines from the list
    for line in error_lines[:10]:
        print(line)
# Example usage:
# response_text = '...'
# if 'ERROR' in response_text:
#     print_error(response_text)

# Note: The actual network request code is omitted as it cannot be executed in this environment.

# Function to download and process files
def download_and_process(nz, ne, entity, prefix):
    filename = create_filename(prefix, nz, ne)
    if os.path.exists(filename):
        return
    base_url = f"http://cdsweb.u-strasbg.fr/cgi-bin/topbase/topbase.sh"
    params = {
        'com': 'dt',
        'ent': entity,
        'nz1': nz,
        'nz2': nz,
        'ne1': ne,
        'ne2': ne,
        'is1': 1,
        'is2': 9,
        'il1': 0,
        'il2': 9,
        'ip1': 0,
        'ip2': 1,
        'lv1': 0,
        'lv2': 0,
        'en1': 0.0,
        'en2': 0.0,
        'so': 's2',
    }
    if entity == 'e':
        params['te'] = 'on'
        params['gi'] = 'on'
    response = requests.get(base_url, params=params)
    
    if response.status_code == 200:
        # Write to file if no error in response
        if 'ERROR' not in response.text:
            with open(filename, 'w') as file:
                file.write(response.text)
            # Remove the first 8 lines of header info
            lines_remove = 8
            with open(filename, 'r') as file:
                lines = file.readlines()
                with open(filename, 'w') as file:
                    file.writelines(lines[lines_remove:])
        else:
            #print(f"*** ERROR ***: nz={nz}, ne={ne}")
            os.remove(filename) # Remove file if there's an error in response
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# Function to wrap the download and process steps
def download_wrapper(nz, ne):
    # Download and process energy level data
    download_and_process(nz, ne, 'e', 'elevel')
    # Download and process photoionization cross section data
    download_and_process(nz, ne, 'p', 'xsectn')

# Determine the total number of tasks to initialize tqdm progress bar
total_tasks = sum(nz for nz in atomic_numbers)

# Using ThreadPoolExecutor to parallelize the downloads
with ThreadPoolExecutor() as executor:
    # Create a list to hold the futures
    futures = []
    for nz in atomic_numbers:
        for ne in range(1, nz + 1):
            # Submit the tasks to the executor
            futures.append(executor.submit(download_wrapper, nz, ne))

    # Iterate over the futures as they complete and update progress bar
    for _ in tqdm(as_completed(futures), total=total_tasks, desc='Downloading', unit='file'):
        pass

# Please note that actual function calls inside download_wrapper are commented out to prevent execution errors in this environment.
# You should implement the download_and_process function with the actual downloading logic in your local environment.


# from tqdm import tqdm

# # ... (the rest of your existing code)

# # Loop through the defined range of atomic numbers and electrons with a progress bar
# for nz in tqdm(atomic_numbers, desc='Atomic numbers', unit='atom'):
#     for ne in tqdm(range(1, nz + 1), desc=f'Electrons for Z={nz}', leave=False, unit='electron'):
#         # Download and process energy level data
#         download_and_process(nz, ne, 'e', 'elevel')

#         # Download and process photoionization cross section data
#         download_and_process(nz, ne, 'p', 'xsectn')

# Convert fetched files to an HDF5 file for Optab
import subprocess
subprocess.call(["../src/convert_topbase_h5"])

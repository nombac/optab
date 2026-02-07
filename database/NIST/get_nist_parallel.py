import csv
import re
import requests
import os
import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed

# Function to process each 'aaa' value
def process_aaa(aaa, output_directory, pbar):
    for bb in range(aaa):
        output_states_path = os.path.join(output_directory, f"nist_{aaa:03d}.{bb:02d}.states")
        output_ionize_path = os.path.join(output_directory, f"nist_{aaa:03d}.{bb:02d}.ionize")

        if os.path.exists(output_states_path) and os.path.exists(output_ionize_path):
            pbar.update(1)
            continue

        url = f"https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=Z%3D{aaa}+{bb}&submit=Retrieve+Data&units=0&format=2&output=0&page_size=15&average_out=1&multiplet_ordered=1&conf_out=on&term_out=on&level_out=on&g_out=on&temp="
        response = requests.get(url)

        if response.status_code == 200:
            output_csv_path = os.path.join(output_directory, f"nist_{aaa:03d}.{bb:02d}.tmp")
            with open(output_csv_path, 'w') as file:
                file.write(response.text)

            with open(output_csv_path, mode='r', newline='') as csv_file, \
                    open(output_states_path, mode='w') as states_file, \
                    open(output_ionize_path, mode='w') as ionize_file:
                csv_reader = csv.reader(csv_file)
                next(csv_reader)  # Skip the header row
                
                pattern = re.compile(r'[="[\]()]+')

                for row in csv_reader:
                    cleaned_row = [pattern.sub('', cell) for cell in row]
                    if "Limit" in cleaned_row:
                        output_row = [cleaned_row[4]]
                        ionize_file.write(" ".join(output_row) + "\n")
                    else:
                        if aaa != 1 or (aaa == 1 and (cleaned_row[0].isdigit() or cleaned_row[0] == "1s")):
                            output_row = [cleaned_row[2], cleaned_row[4], cleaned_row[1], cleaned_row[0]]
                            states_file.write(" ".join(output_row) + "\n")
        else:
            print(f"Error retrieving data for Z={aaa} {bb}: Status code {response.status_code}")

        # Safely update the progress bar
        pbar.update(1)

def main():
    # Set the minimum and maximum values for 'aaa'
    amin = 1
    amax = 92

    # Calculate the total iterations for the progress bar
    total_iterations = sum(range(amin, amax + 1))

    # Set the directory paths for the input CSV and output text files
    output_directory = "levels"

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Initialize the progress bar with the total iterations
    with tqdm.tqdm(total=total_iterations) as pbar:
        # Create a ThreadPoolExecutor to parallelize the 'aaa' loop
        with ThreadPoolExecutor() as executor:
            # Create a future for each 'aaa' in the range
            futures = [executor.submit(process_aaa, aaa, output_directory, pbar) for aaa in range(amin, amax + 1)]
            
            # Wait for all futures to complete
            for future in as_completed(futures):
                future.result()  # Get the result to make sure any exceptions are raised

    # Print completion message
    print("The process has been completed.")

    # Fetch isotope data from NIST (replaces get_nist_atomic.sh / lynx)
    print("Fetching isotope data from NIST...")
    url = "https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&all=all&ascii=ascii2&isotype=all"
    response = requests.get(url)
    if response.status_code == 200:
        text = response.text
        # Extract data records (between first and last "Atomic Number" blocks)
        start = text.find("Atomic Number")
        end = text.find("\n\n", text.rfind("Notes"))
        text = text[start:end]
        text = text.replace("&nbsp;", " ")
        text = text.replace("Atomic Symbol = D", "Atomic Symbol = H")
        text = text.replace("Atomic Symbol = T", "Atomic Symbol = H")
        with open("nist_isotope.txt", "w") as f:
            f.write(text + "\n")
    else:
        print(f"Error fetching isotope data: Status code {response.status_code}")

    # Convert fetched files to an HDF5 file for Optab
    import subprocess
    subprocess.call(["../src/convert_nist_h5"])

if __name__ == "__main__":
    main()

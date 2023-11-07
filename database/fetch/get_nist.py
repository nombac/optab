import csv
import re
import requests
import os
import tqdm

# Set the minimum and maximum values for 'aaa'
amin = 1
amax = 92

# Set the directory paths for the input CSV and output text files
output_directory = "levels"

# Create the output directory if it doesn't exist
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Calculate the total iterations for the progress bar
total_iterations = sum(range(amin, amax + 1))

# Loop through the range of 'aaa' from amin to amax
with tqdm.tqdm(total=total_iterations) as pbar:
    for aaa in range(amin, amax + 1):
        for bb in range(aaa):
            # Generate output file paths using os.path.join
            output_states_path = os.path.join(output_directory, f"nist_{aaa:03d}.{bb:02d}.states")
            output_ionize_path = os.path.join(output_directory, f"nist_{aaa:03d}.{bb:02d}.ionize")

            # Construct the URL for the request
            url = f"https://physics.nist.gov/cgi-bin/ASD/energy1.pl?de=0&spectrum=Z%3D{aaa}+{bb}&submit=Retrieve+Data&units=0&format=2&output=0&page_size=15&average_out=1&multiplet_ordered=1&conf_out=on&term_out=on&level_out=on&g_out=on&temp="
            response = requests.get(url)

            # Check the response status and write to a file if successful
            if response.status_code == 200:
                output_csv_path = os.path.join(output_directory, f"nist_{aaa:03d}.{bb:02d}.tmp")
                with open(output_csv_path, 'w') as file:
                    file.write(response.text)
            else:
                print(f"Error retrieving data for Z={aaa} {bb}: Status code {response.status_code}")

            # Process the CSV file and output the results to text files
            with open(output_csv_path, mode='r', newline='') as csv_file, \
                 open(output_states_path, mode='w') as states_file, \
                 open(output_ionize_path, mode='w') as ionize_file:
                csv_reader = csv.reader(csv_file)
                next(csv_reader)  # Skip the header row

                pattern = re.compile(r'[="[\]()]+')

                for row in csv_reader:
                    # Remove double quotes and equals sign from each line
                    cleaned_row = [pattern.sub('', cell) for cell in row]
                    if "Limit" in cleaned_row:
                        # Write the row containing "Limit" to the ionize file
                        output_row = [cleaned_row[4]]
                        ionize_file.write(" ".join(output_row) + "\n")
                        #print(" ".join(output_row))
                    else:
                        # Keep rows where the first column is only numbers (except for "1s")
                        if aaa != 1 or (aaa == 1 and (cleaned_row[0].isdigit() or cleaned_row[0] == "1s")):
                            # Output in the order of 3rd, 5th, 1st, and 2nd columns
                            output_row = [cleaned_row[2], cleaned_row[4], cleaned_row[1], cleaned_row[0]]
                            states_file.write(" ".join(output_row) + "\n")
                            #print(" ".join(output_row))

            # Update the progress bar
            pbar.set_description(f"Processing aaa={aaa}, bb={bb}")
            pbar.update(1)


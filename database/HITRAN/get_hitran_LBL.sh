#!/bin/bash

# Ensure the directory exists where files will be saved
mkdir -p original

# Array of URLs to download
urls=(
    "https://hitran.org/suppl/LBLstatic/30/30_hit20.zip"
    "https://hitran.org/suppl/LBLstatic/35/35_hit20.zip"
    "https://hitran.org/suppl/LBLstatic/42/42_hit20.zip"
    "https://hitran.org/suppl/LBLstatic/55/55_hit20.zip"
)

# Loop over URLs to download and process each one
for url in "${urls[@]}"; do
    # Extract the file base name (e.g., 30_hit20) from the URL
    base_name=$(basename "$url" .zip)

    # Construct the new .par file name by replacing 'hit20' with 'HITRAN'
    new_par_name="${base_name/hit20/HITRAN}.par"

    # Download the zip file to the original/ directory
    wget -q "$url" -P original/

    # Navigate to the original directory without printing directory stack status
    pushd original/ > /dev/null

    # Unzip the file
    unzip -qo "$base_name.zip"

    # Remove the PDF file
    rm -f "$base_name.pdf"

    # Rename the .par file
    mv "$base_name.par" "$new_par_name"

    # Remove the zip file
    rm "$base_name.zip"

    echo "$base_name done"

    # Return to the previous directory without printing directory stack status
    popd > /dev/null
done

# Inform the user that all files have been processed
echo "All LBL files have been downloaded, extracted, and processed."

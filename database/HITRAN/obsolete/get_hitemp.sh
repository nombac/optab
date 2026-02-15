#!/bin/bash

# Define the URL of the directory containing the bz2 files
BASE_URL="https://hitran.org/hitemp/data/bzip2format/"

# Create a directory where the bz2 files will be downloaded and extracted
mkdir -p original

# Use wget to download all bz2 files from the specified directory into the original directory
wget -r -nd -A '*.bz2' -P original "$BASE_URL"

# Extract all bz2 files in the original directory
find original/ -name '*.bz2' -exec bunzip2 {} \;

import os
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# Function to check if a URL exists
def url_exists(url):
    try:
        response = requests.head(url, allow_redirects=True)
        return response.status_code == 200
    except requests.RequestException as e:
        print(f"Request for {url} caused an exception: {e}")
        return False

# Function to download and save partition function files
def download_partition_function(n):
    fname = f"q{n}.txt"
    #m = f"{n:03d}"
    #oname = f"q{m}.txt"
    oname = f"q{n}.txt"
    URL = f"https://hitran.org/data/Q/{fname}"
    if url_exists(URL):
        response = requests.get(URL)
        if response.status_code == 200:
            with open(os.path.join(dir_path, 'Q', oname), 'wb') as file:
                file.write(response.content)
            return f"Downloaded {oname}"
        else:
            return f"Failed to download {oname}"
    else:
        return f"URL {URL} not found"

dir_path = './'
os.makedirs(os.path.join(dir_path, 'Q'), exist_ok=True)

# Use ThreadPoolExecutor to parallelize downloads
with ThreadPoolExecutor(max_workers=10) as executor:
    # Submit all tasks and store the futures in a list
    futures = {executor.submit(download_partition_function, n): n for n in range(1, 149)}
    
    # Use tqdm to create a progress bar for the futures as they complete
    for future in tqdm(as_completed(futures), total=len(futures), desc="Downloading Q files..."):
        url = futures[future]
        try:
            # Check the result of the future
            result = future.result()
            #print(result)
        except Exception as exc:
            print(f"Partition function {url} generated an exception: {exc}")

# Convert fetched files to an HDF5 file for Optab
import subprocess
subprocess.call(["bash", "../fetch/get_hitran_meta.sh"])


from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import os
import requests

# Function that checks if a URL exists
def url_exists(url):
    try:
        response = requests.head(url, timeout=5)
        return response.status_code == 200
    except requests.RequestException:
        return False

# Function to download a file
def download_file(nz, ne, dir_path):
    number = f"{nz:02d}{ne:02d}"
    ename = f"gf{number}.gam"
    url = f"http://kurucz.harvard.edu/atoms/{number}/{ename}"
    if url_exists(url):
        response = requests.get(url)
        with open(os.path.join(dir_path, ename), 'wb') as f:
            f.write(response.content)
        return f"{ename} is fetched..."
    else:
        return f"{ename} does not exist..."

# Set directory and range for nz
directory = './atoms/'
nz_min = 1
nz_max = 92

# Ensure the directory exists
if not os.path.exists(directory):
    os.makedirs(directory)

# Prepare the tasks
tasks = []
for nz in range(nz_min, nz_max + 1):
    ne_max = nz - 1
    for ne in range(ne_max + 1):
        tasks.append((nz, ne, directory))

# Execute the tasks in parallel and display a progress bar
with ThreadPoolExecutor() as executor:
    futures = [executor.submit(download_file, nz, ne, directory) for nz, ne, directory in tasks]
    results = []
    for future in tqdm(as_completed(futures), total=len(futures), desc="Downloading", unit="file"):
        results.append(future.result())

# Results contain messages about fetched files or errors
for result in results:
    print(result)


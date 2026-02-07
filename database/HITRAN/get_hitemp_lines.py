"""
Download HITEMP line-by-line data for all available molecules.
Uses Chrome user profile for HITRAN authentication.

HITEMP molecules:
  01: H2O  (multi-part zip, HITEMP-2010)
  02: CO2  (bz2, HITEMP2024)
  04: N2O  (bz2, HITEMP2019)
  05: CO   (bz2, HITEMP2019)
  06: CH4  (bz2, HITEMP2020)
  08: NO   (bz2, HITEMP2019)
  10: NO2  (bz2, HITEMP2019)
  13: OH   (bz2, HITEMP2020)

Usage:
  python3 get_hitemp_lines.py "/path/to/Chrome/Profile"
"""

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager
import requests
import os
import bz2
import zipfile
import shutil
import argparse

parser = argparse.ArgumentParser(description='Download HITEMP data.')
parser.add_argument('user_profile_dir', help='The path to the Chrome user profile directory')
args = parser.parse_args()

USER_PROFILE_DIR = args.user_profile_dir

# HITEMP molecules available as single bz2 files
bz2_files = [
    ('02_HITEMP2024.par.bz2', '02_HITEMP.par'),
    ('04_HITEMP2019.par.bz2', '04_HITEMP.par'),
    ('05_HITEMP2019.par.bz2', '05_HITEMP.par'),
    ('06_HITEMP2020.par.bz2', '06_HITEMP.par'),
    ('08_HITEMP2019.par.bz2', '08_HITEMP.par'),
    ('10_HITEMP2019.par.bz2', '10_HITEMP.par'),
    ('13_HITEMP2020.par.bz2', '13_HITEMP.par'),
]

BZ2_BASE_URL = 'https://hitran.org/files/HITEMP/bzip2format/'
H2O_BASE_URL = 'https://hitran.org/files/HITEMP/HITEMP-2010/H2O_line_list/'

os.makedirs('original', exist_ok=True)


def get_authenticated_session(driver):
    """Extract cookies from Selenium and create an authenticated requests session."""
    session = requests.Session()
    for cookie in driver.get_cookies():
        session.cookies.set(cookie['name'], cookie['value'])
    return session


def download_streaming(session, url, filepath):
    """Stream-download a file with progress indication."""
    response = session.get(url, stream=True)
    if response.status_code != 200:
        print(f'  Failed with status code: {response.status_code}')
        return False
    total = int(response.headers.get('content-length', 0))
    downloaded = 0
    with open(filepath, 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)
            downloaded += len(chunk)
            if total > 0:
                pct = downloaded * 100 // total
                print(f'\r  {pct}% ({downloaded // (1024*1024)} / {total // (1024*1024)} MB)',
                      end='', flush=True)
    print()
    return True


# Set up Chrome options
chrome_options = Options()
chrome_options.add_argument(f"user-data-dir={USER_PROFILE_DIR}")

# ======================================================================
# Phase 1: Download bz2-compressed files (CO2, N2O, CO, CH4, NO, NO2, OH)
# ======================================================================

# Check if any bz2 files need downloading
bz2_to_download = [(bz2_name, par_name) for bz2_name, par_name in bz2_files
                    if not os.path.exists(os.path.join('original', par_name))]

if bz2_to_download:
    print("=== Downloading bz2-compressed HITEMP files ===")
    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=chrome_options)

    try:
        # Navigate to authenticate
        driver.get(BZ2_BASE_URL)
        WebDriverWait(driver, 300).until(
            EC.presence_of_element_located((By.PARTIAL_LINK_TEXT, ".par.bz2"))
        )
        print("Authentication successful.")

        session = get_authenticated_session(driver)

        for bz2_name, par_name in bz2_to_download:
            filepath = os.path.join('original', par_name)
            bz2_path = os.path.join('original', bz2_name)

            print(f'Downloading {bz2_name}...')
            if download_streaming(session, BZ2_BASE_URL + bz2_name, bz2_path):
                print(f'  Decompressing...')
                with bz2.BZ2File(bz2_path, 'rb') as bz2_in:
                    with open(filepath, 'wb') as par_out:
                        shutil.copyfileobj(bz2_in, par_out)
                os.remove(bz2_path)
                print(f'  Saved to {filepath}.')
    finally:
        driver.quit()
else:
    for _, par_name in bz2_files:
        print(f'Skipping original/{par_name} (already exists).')

# ======================================================================
# Phase 2: Download H2O multi-part zip files
# ======================================================================

filepath_h2o = os.path.join('original', '01_HITEMP.par')
if os.path.exists(filepath_h2o):
    print(f'Skipping {filepath_h2o} (already exists).')
else:
    print("\n=== Downloading H2O HITEMP data (multi-part zip) ===")
    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=chrome_options)

    try:
        driver.get(H2O_BASE_URL)
        WebDriverWait(driver, 300).until(
            EC.presence_of_element_located((By.PARTIAL_LINK_TEXT, ".zip"))
        )
        print("Authentication successful.")

        session = get_authenticated_session(driver)

        # Collect all zip links and sort by wavenumber range
        links = driver.find_elements(By.TAG_NAME, "a")
        zip_urls = sorted([
            link.get_attribute("href")
            for link in links
            if link.get_attribute("href") and link.get_attribute("href").endswith(".zip")
        ])
        print(f'Found {len(zip_urls)} zip files.')

        # Download, extract and concatenate into a single .par file
        with open(filepath_h2o, 'wb') as concat_file:
            for i, url in enumerate(zip_urls):
                zip_name = os.path.basename(url)
                zip_path = os.path.join('original', zip_name)
                print(f'  [{i+1}/{len(zip_urls)}] {zip_name}')

                if download_streaming(session, url, zip_path):
                    with zipfile.ZipFile(zip_path, 'r') as zf:
                        for name in sorted(zf.namelist()):
                            if name.endswith('.par'):
                                with zf.open(name) as par_file:
                                    shutil.copyfileobj(par_file, concat_file)
                    os.remove(zip_path)

        print(f'Saved to {filepath_h2o}.')
    finally:
        driver.quit()

print("\nAll HITEMP files downloaded.")

"""
NOTE
If this code does not work, download the files manually as follows (e.g. H2O; repeat this procedure for other species):

1. Goto [Line-by-Line Search](https://hitran.org/lbl/).
1. "Select Molecules" &rarr; check 1. H2O
1. "Select Isotopologues" &rarr; check all isotopologues
1. "Select Wavenumber / Wavelength Range" &rarr; leave blank for &nu;<sub>max</sub>
1. "Select or Create Output Format" &rarr; .par (160 chars)
1. "Start Data Search> Search Results" &rarr; download the "Output transitions data (160-character `.par` format)" as `original/01_HITRAN.par`. Here, `"01"` is the two digits [molecule ID](https://hitran.org/docs/molec-meta/) of H<sub>2</sub>O.
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
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Download HITRAN data.')
parser.add_argument('user_profile_dir', help='The path to the Chrome user profile directory')

# Parse arguments
args = parser.parse_args()

USER_PROFILE_DIR = args.user_profile_dir

# List of (URL, filename) tuples
urls_and_files = [
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=1%2C2%2C3%2C4%2C5%2C6%2C129&vib_bands=&numin=0&numax=', '01_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=7%2C8%2C9%2C10%2C11%2C12%2C13%2C14%2C15%2C121&vib_bands=&numin=0&numax=', '02_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=16%2C17%2C18%2C19%2C20&vib_bands=&numin=0&numax=', '03_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=21%2C22%2C23%2C24%2C25&vib_bands=&numin=0&numax=', '04_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=26%2C27%2C28%2C29%2C30%2C31&vib_bands=&numin=0&numax=', '05_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=32%2C33%2C34%2C35&vib_bands=&numin=0&numax=', '06_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=36%2C37%2C38&vib_bands=&numin=0&numax=', '07_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=39%2C40%2C41&vib_bands=&numin=0&numax=', '08_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=42%2C43%2C137%2C138&vib_bands=&numin=0&numax=', '09_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=44%2C130&vib_bands=&numin=0&numax=', '10_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=45%2C46&vib_bands=&numin=0&numax=', '11_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=47%2C117&vib_bands=&numin=0&numax=', '12_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=48%2C49%2C50&vib_bands=&numin=0&numax=', '13_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=51%2C110&vib_bands=&numin=0&numax=', '14_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=52%2C53%2C107%2C108&vib_bands=&numin=0&numax=', '15_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=54%2C55%2C111%2C112&vib_bands=&numin=0&numax=', '16_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=56%2C113&vib_bands=&numin=0&numax=', '17_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=57%2C58&vib_bands=&numin=0&numax=', '18_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=59%2C60%2C61%2C62%2C63%2C135&vib_bands=&numin=0&numax=', '19_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=64%2C65%2C66&vib_bands=&numin=0&numax=', '20_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=67%2C68&vib_bands=&numin=0&numax=', '21_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=69%2C118&vib_bands=&numin=0&numax=', '22_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=70%2C71%2C72&vib_bands=&numin=0&numax=', '23_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=73%2C74&vib_bands=&numin=0&numax=', '24_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=75&vib_bands=&numin=0&numax=', '25_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=76%2C77%2C105&vib_bands=&numin=0&numax=', '26_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=78%2C106&vib_bands=&numin=0&numax=', '27_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=79&vib_bands=&numin=0&numax=', '28_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=80%2C119&vib_bands=&numin=0&numax=', '29_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=81%2C82%2C83&vib_bands=&numin=0&numax=', '31_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=84&vib_bands=&numin=0&numax=', '32_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=85&vib_bands=&numin=0&numax=', '33_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=86&vib_bands=&numin=0&numax=', '34_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=87&vib_bands=&numin=0&numax=', '36_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=88%2C89&vib_bands=&numin=0&numax=', '37_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=90%2C91&vib_bands=&numin=0&numax=', '38_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=92&vib_bands=&numin=0&numax=', '39_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=93%2C94&vib_bands=&numin=0&numax=', '40_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=95&vib_bands=&numin=0&numax=', '41_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=116&vib_bands=&numin=0&numax=', '43_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=109&vib_bands=&numin=0&numax=', '44_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=103%2C115&vib_bands=&numin=0&numax=', '45_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=97%2C98%2C99%2C100&vib_bands=&numin=0&numax=', '46_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=114&vib_bands=&numin=0&numax=', '47_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=123&vib_bands=&numin=0&numax=', '48_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=124%2C125&vib_bands=&numin=0&numax=', '49_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=146%2C147%2C148&vib_bands=&numin=0&numax=', '50_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=144&vib_bands=&numin=0&numax=', '51_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=139%2C140%2C141%2C142%2C143&vib_bands=&numin=0&numax=', '52_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=131%2C132%2C133%2C134&vib_bands=&numin=0&numax=', '53_HITRAN.par'),
    ('https://hitran.org/lbl/5?output_format_id=1&iso_ids_list=145&vib_bands=&numin=0&numax=', '54_HITRAN.par')
    # Add more tuples here for each URL and target file name
]

# Set up Chrome options
chrome_options = Options()
chrome_options.add_argument(f"user-data-dir={USER_PROFILE_DIR}")

# Ensure the directory exists where files will be saved
os.makedirs('original', exist_ok=True)

# Iterate over the URL and filename pairs
for URL, parfile in urls_and_files:
    # Set up the WebDriver with a new Service object each time
    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=chrome_options)
    # Open the URL
    driver.get(URL)

    # Wait for the file link to be generated
    try:
        # Wait for the .par file link to appear (adjust the timeout as needed)
        file_link = WebDriverWait(driver, 180).until(
            EC.presence_of_element_located((By.PARTIAL_LINK_TEXT, ".par"))
        )
        file_url = file_link.get_attribute('href')

        # Output the URL of the .par file
        print(f'File URL found: {file_url}')

        # Download the file using requests
        response = requests.get(file_url)
        if response.status_code == 200:
            filepath = os.path.join('original', parfile)
            with open(filepath, 'wb') as f:
                f.write(response.content)
            print(f'File downloaded successfully to {filepath}.')
        else:
            print('Failed to download the file.')

    finally:
        driver.quit()

import subprocess        
subprocess.run(['bash', 'get_hitran_LBL.sh'])

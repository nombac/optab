from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import os
import time
import zipfile
import shutil

# 現在のワーキングディレクトリを取得し、相対パスでダウンロード先を指定
current_working_dir = os.getcwd()  # 実行しているスクリプトのワーキングディレクトリを取得
download_dir = os.path.join(current_working_dir, "original")  # "./original" を絶対パス化
os.makedirs(download_dir, exist_ok=True)  # ダウンロード先ディレクトリを作成（存在しない場合）

# Chromeのオプション設定
options = webdriver.ChromeOptions()
prefs = {
    "profile.default_content_settings.popups": 0,  # ポップアップを無効化
    "download.default_directory": download_dir,   # ダウンロード先ディレクトリ
    "directory_upgrade": True,                    # ディレクトリのアップグレードを許可
    "download.prompt_for_download": False,        # ダウンロードプロンプトを無効化
    "safebrowsing.enabled": True                  # セーフブラウジングを有効化
}
options.add_experimental_option("prefs", prefs)

# ChromeDriverの起動
driver = webdriver.Chrome(options=options)

def wait_for_downloads(directory, timeout=1800):
    """指定されたディレクトリ内のダウンロードが完了するまで待機する"""
    start_time = time.time()
    while True:
        # ダウンロード中のファイルが存在するか確認
        downloading_files = [f for f in os.listdir(directory) if f.endswith(".crdownload")]
        if not downloading_files:  # ダウンロード中のファイルがなければ終了
            break
        if time.time() - start_time > timeout:  # タイムアウト処理
            raise TimeoutError("Downloads did not complete within the given time.")
        time.sleep(2)  # 2秒待機して再確認

def decompress_zip(directory):
    """ディレクトリ内の.zipファイルを解凍し、元ファイルを削除する"""
    for filename in os.listdir(directory):
        if filename.endswith(".zip"):
            filepath = os.path.join(directory, filename)
            decompressed_dir = os.path.join(directory, filename[:-4])  # .zipを取り除いた名前のディレクトリ
            print(f"Decompressing: {filepath}")
            with zipfile.ZipFile(filepath, 'r') as zip_ref:
                zip_ref.extractall(decompressed_dir)  # 解凍先ディレクトリを指定
            os.remove(filepath)  # 解凍後に元ファイルを削除
            print(f"Deleted: {filepath}")

def concatenate_and_cleanup(directory):
    """展開された.parファイルを結合し、元ファイルを削除する"""
    concat_file_path = os.path.join(directory, "01_HITEMP2010.par")
    with open(concat_file_path, "wb") as concat_file:
        subdirs = [d for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d)) and d.startswith("01_")]
        sorted_subdirs = sorted(subdirs, key=lambda x: int(x.split("_")[1].split("-")[0]))
        for subdir in sorted_subdirs:
            subdir_path = os.path.join(directory, subdir)
            par_files = [f for f in os.listdir(subdir_path) if f.endswith(".par")]
            for par_file in par_files:
                par_file_path = os.path.join(subdir_path, par_file)
                print(f"Adding {par_file_path} to {concat_file_path}")
                with open(par_file_path, "rb") as pf:
                    shutil.copyfileobj(pf, concat_file)
            # サブディレクトリとその中のファイルを削除
            shutil.rmtree(subdir_path)
            print(f"Deleted directory: {subdir_path}")

try:
    # HITRANページにアクセス
    url = "https://hitran.org/files/HITEMP/HITEMP-2010/H2O_line_list/"
    driver.get(url)
    print("Please log in manually in the browser...")

    # ログイン完了を待機
    WebDriverWait(driver, 300).until(
        EC.presence_of_element_located((By.XPATH, "//a[contains(@href, '.zip')]"))
    )
    print("Login detected. Proceeding with downloads...")

    # ページ内のリンクを取得
    links = driver.find_elements(By.TAG_NAME, "a")

    # .zip で終わるリンクをクリックしてダウンロード
    for link in links:
        href = link.get_attribute("href")
        if href and href.endswith(".zip"):
            print(f"Downloading: {href}")
            link.click()
            time.sleep(5)  # ダウンロードが開始されるのを待つ

    # ダウンロードが完了するまで待機
    print("Waiting for downloads to complete...")
    wait_for_downloads(download_dir)
    print("All downloads complete.")

finally:
    # ブラウザを閉じる
    driver.quit()

# ダウンロードした.zipファイルを解凍し、元ファイルを削除
print("Decompressing downloaded files...")
decompress_zip(download_dir)

# 展開された.parファイルを結合し、元ファイルを削除
print("Concatenating .par files...")
concatenate_and_cleanup(download_dir)
print("All files concatenated and original directories removed.")

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import os
import time
import bz2
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

def decompress_bz2(directory):
    """ディレクトリ内の.bz2ファイルを解凍し、元ファイルを削除する"""
    for filename in os.listdir(directory):
        if filename.endswith(".bz2"):
            filepath = os.path.join(directory, filename)
            decompressed_path = os.path.join(directory, filename[:-4])  # .bz2を取り除いた名前
            print(f"Decompressing: {filepath}")
            with bz2.BZ2File(filepath, 'rb') as bz2_file:
                with open(decompressed_path, 'wb') as decompressed_file:
                    shutil.copyfileobj(bz2_file, decompressed_file)
            os.remove(filepath)  # 解凍後に元ファイルを削除
            print(f"Deleted: {filepath}")

try:
    # HITRANページにアクセス
    url = "https://hitran.org/files/HITEMP/bzip2format/"
    driver.get(url)
    print("Please log in manually in the browser...")

    # ログイン完了を待機
    WebDriverWait(driver, 300).until(
        EC.presence_of_element_located((By.XPATH, "//a[contains(@href, '.par.bz2')]"))
    )
    print("Login detected. Proceeding with downloads...")

    # ページ内のリンクを取得
    links = driver.find_elements(By.TAG_NAME, "a")

    # .par.bz2 で終わるリンクをクリックしてダウンロード
    for link in links:
        href = link.get_attribute("href")
        if href and href.endswith(".par.bz2"):
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

# ダウンロードした.bz2ファイルを解凍し、元ファイルを削除
print("Decompressing downloaded files...")
decompress_bz2(download_dir)
print("All files decompressed and original .bz2 files deleted.")

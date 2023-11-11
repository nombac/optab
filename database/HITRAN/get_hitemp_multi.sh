#!/bin/bash

download_and_process() {
  id=$1
  directory=$2
  name="HITEMP2010"

  # データをダウンロード
  wget -r -l1 -H -t1 -nd -N -np -A.zip -erobots=off https://hitran.org/hitemp/data/${directory}

  # カレントディレクトリ内の全ての .zip ファイルに対してループ処理
  for zipfile in *.zip; do
    # 拡張子を除いたファイル名を取得
    base=${zipfile%.zip}
    
    # 一時ディレクトリを作成
    tmpdir=$(mktemp -d)
    
    # .zip ファイルを一時ディレクトリに解凍
    unzip "$zipfile" -d "$tmpdir"
    
    # 一時ディレクトリ内の .par ファイルを探し、新しいファイル名でカレントディレクトリに移動
    find "$tmpdir" -name '*.par' -exec mv {} "${base}.par" \;
    
    # 一時ディレクトリを削除
    rm -r "$tmpdir" "$zipfile"
  done

  # .par ファイルを結合
  cat ${id}_?????-?????_${name}.par > ${id}_${name}.par
  # 不要な .par ファイルを削除
  rm -r ${id}_?????-?????_${name}.par
}

mkdir -p original
cd original/

# H2O のデータを処理
download_and_process "01" "HITEMP-2010/H2O_line_list/"

# CO2 のデータを処理
download_and_process "02" "HITEMP-2010/CO2_line_list/"

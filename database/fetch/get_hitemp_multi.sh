#!/bin/bash

name="HITEMP2010"
id="01"
directory="HITEMP-2010/H2O_line_list/"

cd original

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

cat ${id}_?????-?????_${name}.par > ${id}_${name}.par
rm -r ${id}_?????-?????_${name}.par

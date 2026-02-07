#!/bin/bash

mkdir -p decomp

for file in original/*_HITEMP.par
do
    echo ""
    echo "*** $file ****"
    echo "preprocessing..."
    ../src/preproc_hitran $file
    echo "converting..."
    ../src/convert_lines_h5
    echo "done"
done

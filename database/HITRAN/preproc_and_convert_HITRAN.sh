#!/bin/bash

for file in original/*.par
do
    echo $file...
    echo preprocessing...
    ../src/preproc_hitran $file
    echo converting...
    ../src/convert_lines_h5    
done

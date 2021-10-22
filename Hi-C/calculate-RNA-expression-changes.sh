#!/bin/bash

for dir in ESC 2CLC
do

for file in AtoB.bed BtoA.bed
do
    file2=${file##*/}
    ./cal_rpkm.sh ${dir} ${file} repeat_expression50/${dir}-${file2/.bed/}.rpkm
done
done

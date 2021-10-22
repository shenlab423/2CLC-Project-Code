#!/bin/bash

mappingfile=$1
outfile=$2

/nethome/yezhang_zhu/miniconda3/bin/samtools view -@ 10 -H $mappingfile > $outfile
/nethome/yezhang_zhu/miniconda3/bin/samtools view -@ 10 $mappingfile | LC_ALL=C fgrep -w "NH:i:1" >> $outfile

#/nethome/yezhang_zhu/miniconda3/bin/samtools view -@ 10 $mappingfile | grep -w "NH:i:1" >> $outfile

/nethome/yezhang_zhu/miniconda3/bin/samtools view -@ 10 -bh $outfile > ${outfile}.bam

rm $outfile

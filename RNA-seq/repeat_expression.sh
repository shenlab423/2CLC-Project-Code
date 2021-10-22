for dir in `ls *bam`
do
for file in mm9_repeats_simple.bed
do
    ./cal_rpkm.sh ${dir} ${file} repeat_expression/${dir/.unique.bam/}-${file/.bed/}.rpkm 
done
done


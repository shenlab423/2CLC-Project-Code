neg=$1
pos=$2
name=$3
neg_sam=$4
pos_sam=$5

cat $neg $pos > $name_peaks.broadPeak
bedtools sort -i $name_peaks.broadPeak > $name_peaks.sort.broadPeak
bedtools merge -i $name_peaks.sort.broadPeak > $name_peaks.merge.broadPeak
bash cal_rpkm.sh $neg_sam $name_peaks.merge.broadPeak
bash cal_rpkm.sh $pos_sam $name_peaks.merge.broadPeak
paste ${neg_sam##*/}.peaks.IP.rpkm.bdg ${pos_sam##*/}.peaks.IP.rpkm.bdg > ${name}_ESC_2C_rpkm.bdg

awk -v OFS="\t" '{if($4>1||$8>1) print $1,$2,$3,$4,$8}' ${name}_ESC_2C_rpkm.bdg > ${name}_ESC_2C_rpkm.filtered.bdg
awk -v OFS="\t" '{if($4>2*$5) print $1,$2,$3,$4,$5}' ${name}_ESC_2C_rpkm.filtered.bdg > ${name}_2C_DOWN_rpkm.bdg
awk -v OFS="\t" '{if($5>2*$4) print $1,$2,$3,$4,$5}' ${name}_ESC_2C_rpkm.filtered.bdg > ${name}_2C_UP_rpkm.bdg
bedtools intersect -v -a ${name}_ESC_2C_rpkm.filtered.bdg -b ${name}_2C_UP_rpkm.bdg >1.bdg
bedtools intersect -v -a 1.bdg -b ${name}_2C_DOWN_rpkm.bdg > common_peaks.bdg
cat common_peaks.bdg ${name}_2C_DOWN_rpkm.bdg > ${name}_ESC_peaks.bed
cat common_peaks.bdg ${name}_2C_UP_rpkm.bdg > ${name}_2C_peaks.bed
mkdir FinalPeak
mv ${name}_ESC_peaks.bed FinalPeak
mv ${name}_2C_peaks.bed FinalPeak


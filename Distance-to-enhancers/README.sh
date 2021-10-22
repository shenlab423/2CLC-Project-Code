# get Up,Down,Non-regulated promoters using GetUpDownNon-regulated-promoters.R

bedtools sort -i Up_promoter.bed > Up_promoter.sort.bed
bedtools closest -D ref -t first -a Up_promoter.sort.bed -b OSN.sort.bed > Up_promoter_nearest_OSN.txt

bedtools sort -i Down_promoter.bed > Down_promoter.sort.bed
bedtools closest -D ref -t first -a Down_promoter.sort.bed -b OSN.sort.bed > Down_promoter_nearest_OSN.txt

bedtools sort -i Nochange_promoter.bed > Nochange_promoter.sort.bed
bedtools closest -D ref -t first -a Nochange_promoter.sort.bed -b OSN.sort.bed > Nochange_promoter_nearest_OSN.txt

# calculate distance of ESC enhancers to three type of promoters using distance-Fig2E.R


bedtools sort -i H3K27ac_2C_UP_rpkm.intergenic.bed > H3K27ac_2C_UP_rpkm.intergenic.sort.bed
bedtools closest -D ref -t first -a Up_promoter.sort.bed -b H3K27ac_2C_UP_rpkm.intergenic.sort.bed > Up_promoter_nearest_2Cenhancer.txt
bedtools closest -D ref -t first -a Down_promoter.sort.bed -b H3K27ac_2C_UP_rpkm.intergenic.sort.bed > Down_promoter_nearest_2Cenhancer.txt
bedtools closest -D ref -t first -a Nochange_promoter.sort.bed -b H3K27ac_2C_UP_rpkm.intergenic.sort.bed > Nochange_promoter_nearest_2Cenhancer.txt

# calculate distance of putative 2C enhancers to three type of promoters using distance-Fig3E.R

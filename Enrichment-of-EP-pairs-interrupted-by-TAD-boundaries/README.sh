# Get Top10% up-regulated, Top30% up-regulated, Up-regulated, Down-regulated, Non-regulated promoters using GetgenesPromoter.R

bedtools sort -i Up30_promoter_enrich.bed > Up30_promoter.sort.bed
bedtools sort -i Up10_promoter_enrich.bed > Up10_promoter.sort.bed
bedtools sort -i H3K27ac_2C_UP_rpkm_FC4.intergenic.bed > H3K27ac_2C_UP_rpkm_FC4.intergenic.sort.bed
bedtools sort -i Up_promoter_enrich.bed > Up_promoter.sort.bed
bedtools sort -i Down_promoter_enrich.bed > Down_promoter.sort.bed
bedtools sort -i Nochange_promoter_enrich.bed > Nochange_promoter.sort.bed

bedtools closest -D ref -t first -a Up_promoter.sort.bed -b H3K27ac_2C_UP_rpkm_FC4.intergenic.sort.bed > Up.close.2Cenhancer.len
bedtools closest -D ref -t first -a Down_promoter.sort.bed -b H3K27ac_2C_UP_rpkm_FC4.intergenic.sort.bed > Down.close.2Cenhancer.len
bedtools closest -D ref -t first -a Nochange_promoter.sort.bed -b H3K27ac_2C_UP_rpkm_FC4.intergenic.sort.bed > Nochange.close.2Cenhancer.len
bedtools closest -D ref -t first -a Up30_promoter.sort.bed -b H3K27ac_2C_UP_rpkm_FC4.intergenic.sort.bed > Up30.close.2Cenhancer.len
bedtools closest -D ref -t first -a Up10_promoter.sort.bed -b H3K27ac_2C_UP_rpkm_FC4.intergenic.sort.bed > Up10.close.2Cenhancer.len

# get EP pairs and Control pairs to calculate enrichment
Rscript control.R

rm enrich.txt

file=mESC_Dixon2012-raw_TADs.boundary.bed

echo -e $(bedtools intersect -u -wa -a Down.close.Dixon.control.txt -b $file | wc -l)"\t"$(wc -l Down.close.Dixon.control.txt)"\t"$(bedtools intersect -u -wa -a Down.close.Dixon.txt -b $file | wc -l)"\t"$(wc -l Down.close.Dixon.txt) >> enrich.txt


echo -e $(bedtools intersect -u -wa -a Up.close.Dixon.control.txt -b $file | wc -l)"\t"$(wc -l Up.close.Dixon.control.txt)"\t"""$(bedtools intersect -u -wa -a Up.close.Dixon.txt -b $file | wc -l)"\t"$(wc -l Up.close.Dixon.txt) >> enrich.txt


echo -e $(bedtools intersect -u -wa -a Nochange.close.Dixon.control.txt -b $file | wc -l)"\t"$(wc -l Nochange.close.Dixon.control.txt)"\t"$(bedtools intersect -u -wa -a Nochange.close.Dixon.txt -b $file | wc -l)"\t"$(wc -l Nochange.close.Dixon.txt) >> enrich.txt


echo -e $(bedtools intersect -u -wa -a Top.close.Dixon.control.txt -b $file | wc -l)"\t"$(wc -l Top.close.Dixon.control.txt)"\t"$(bedtools intersect -u -wa -a Top.close.Dixon.txt -b $file | wc -l)"\t"$(wc -l Top.close.Dixon.txt) >> enrich.txt

echo -e $(bedtools intersect -u -wa -a Top30.close.Dixon.control.txt -b $file | wc -l)"\t"$(wc -l Top30.close.Dixon.control.txt)"\t"$(bedtools intersect -u -wa -a Top30.close.Dixon.txt -b $file | wc -l)"\t"$(wc -l Top30.close.Dixon.txt) >> enrich.txt

# generate Fig 3G using Fig3G.R

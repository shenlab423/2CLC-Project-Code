bedtools shuffle -i OSN.bed -chrom -g ~/tools/ucsc_tools/mm9.chrom.sizes > OSN.shuffle.bed


a=$(bedtools intersect -wa -a OSN.bed -b AtoB.bed | wc -l)
b=$(bedtools intersect -wa -a OSN.shuffle.bed -b AtoB.bed | wc -l)
c=$(bedtools intersect -wa -a OSN.bed -b BtoA.bed | wc -l)
d=$(bedtools intersect -wa -a OSN.shuffle.bed -b BtoA.bed | wc -l)

echo -e "$a\t$b\t$c\t$d"

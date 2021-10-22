#!/bin/bash

GTF=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf
GENOME_FA=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome.fa
GENOME=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome

for file in `ls *R1*.fastq`
do 
{
/nethome/yezhang_zhu/tools/mycode/trim-c $file ${file%%_*}_1.fastq 0 50
} &
done
wait

mkdir backup_fastq
mv *R1*.fastq backup_fastq

for file in `ls *_1.fastq`
do
outdir=${file%_*}_tophat
/nethome/yezhang_zhu/tools/tophat-2.1.1.Linux_x86_64/tophat -o ${outdir} --no-coverage-search --no-novel-juncs -p 16 -G ${GTF} ${GENOME} ${file}
done

for filedir in `ls -F | grep 'tophat/$'`
do
{
#echo $filedir
bash ~/tools/mycode/filter_unique_reads_from_tophat.sh ${filedir}accepted_hits.bam ${filedir}accepted_hits_unique
}&
done
wait

for filedir in `ls -F | grep 'tophat/$'`
do
outdir=${filedir%%/}_cuffquant
/nethome/yezhang_zhu/tools/cufflinks-2.2.1.Linux_x86_64/cuffquant -o ${outdir} -p 16 -u $GTF ${filedir}accepted_hits_unique.bam
done

files2=""
label=""
for file in `ls -F | grep '_cuffquant/$'`
do
label=${label}${file%%_*}","
files2=${files2}${file}"abundances.cxb "
done

label=${label%,*}

/nethome/yezhang_zhu/tools/cufflinks-2.2.1.Linux_x86_64/cuffnorm -p 12 --labels ${label} --library-norm-method=classic-fpkm -o cuffnorm_table $GTF ${files2}

mkdir tophat cuffquant
mv *_tophat/ tophat
mv *_cuffquant/ cuffquant

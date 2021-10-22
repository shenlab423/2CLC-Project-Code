for file in *fastq
do
{
trim_galore $file
}&
done
wait

for file in *trimmed.fq; do bowtie2 -p 20 --local --very-sensitive -I 10 -x /nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome -S ${file%%_*}.sam -U $file &> ${file%%_*}.log ; done

for file in *sam
do
samtools sort -@ 10 $file > ${file/.sam/.bam}
samtools rmdup -s ${file/.sam/.bam} ${file/.sam/.rmdup.bam}
done

for file in *rmdup.bam; do macs2 callpeak -f BAM -q 0.05 --nomodel --extsize 120 -B --SPMR -g mm --outdir ${file/.bam/} -n ${file/.bam/} -t $file; done

for dir in `ls -d *`; do file=$dir/${dir}_treat_pileup.bdg; wigToBigWig -clip $file ~/tools/ucsc_tools/mm9.chrom.sizes ${file/.bdg/.bw}; done

for file in *rmdup.bam; do java -jar ~/tools/gem/gem.jar --d /nethome/yezhang_zhu/tools/gem/Read_Distribution_default.txt --g /nethome/yezhang_zhu/tools/gem/mm9.chrom.sizes --genome /nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa --expt $file --f BAM --out ${file/.bam/.GEM} --k_min 6 --k_max 13 --t 20; done


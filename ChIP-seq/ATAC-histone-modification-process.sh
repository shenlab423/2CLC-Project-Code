GENOME=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome


for file in `ls *R1*.fastq.gz`
do
{
trim-gz-c $file ${file%%_*}_1.fastq.gz 0 50
}&
done
wait

for file in `ls *_1.fastq.gz`
do
bowtie2 -p 6 -x $GENOME -U $file -S ${file/_1.fastq.gz/}.sam 2> ${file/_1.fastq.gz/}.align.log
done


for file in `ls *.sam`
do
file2=${file#*/}
file3=${file2/.sam/}
macs2 callpeak -f SAM -q 0.05 --nomodel --nolambda --broad --extsize 300 -B --SPMR -g mm --outdir $file3 -n $file3 -t $file
done

for dir in `ls -d *`
do
file=$dir/${dir}_treat_pileup.bdg
wigToBigWig -clip $file ~/tools/ucsc_tools/mm9.chrom.sizes ${file/.bdg/.bw}
done

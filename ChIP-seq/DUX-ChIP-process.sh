trim_galore Unsorted-plusDox-HAChIP-12hr-rep1_R1.fastq
bowtie2 -p 20 --local --very-sensitive -I 10 -x /nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome -S Unsorted-plusDox-HAChIP-12hr-rep1.sam -U Unsorted-plusDox-HAChIP-12hr-rep1_R1_trimmed.fq

grep -v "XS:i:" Unsorted-plusDox-HAChIP-12hr-rep1.sam > Unsorted-plusDox-HAChIP-12hr-rep1.unique.sam

macs2 callpeak -f SAM -q 0.05 --nomodel --extsize 120 -B --SPMR -g mm --outdir Unsorted-plusDox-HAChIP-12hr-rep1.unique -n Unsorted-plusDox-HAChIP-12hr-rep1.unique -t Unsorted-plusDox-HAChIP-12hr-rep1.unique.sam
cd Unsorted-plusDox-HAChIP-12hr-rep1.unique
wigToBigWig -clip Unsorted-plusDox-HAChIP-12hr-rep1.unique_treat_pileup.bdg ~/tools/ucsc_tools/mm9.chrom.sizes Unsorted-plusDox-HAChIP-12hr-rep1.unique_treat_pileup.bw

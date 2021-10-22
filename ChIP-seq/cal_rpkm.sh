#!/bin/bash

if [ $# -lt 2 ];then
echo "Need 2 parameters! <IP_sample> <peaklist>"
exit
fi

EXTEND=300

IP_sample=$1
IP_mappable_read_count=$(samtools view -F 0x0004 ${IP_sample}.sort.bam | wc -l)


bamToBed -i ${IP_sample}.sort.bam | extend_single2 ~/tools/ucsc_tools/mm9.chrom.sizes $EXTEND|sort -k1,1 -k2,2n > ${IP_sample}.bed

peakfile=$2

cut -f 1,2,3 ${peakfile} | coverageBed -b ${IP_sample}.bed -a - | awk -v OFS='\t' -v SIZE=$IP_mappable_read_count '{print $1,$2,$3,($4*1000000/SIZE)*1000/($3-$2)}' > ${IP_sample}.peaks.IP.rpkm.bdg


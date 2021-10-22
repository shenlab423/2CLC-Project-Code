if [ $# -lt 3 ];then
echo "Need 2 parameters! <IP_sample> <peaklist> <outfile>"
exit
fi

IP_sample=$1
outfile=$3
IP_mappable_read_count=$(samtools view -F 0x0004 $IP_sample | wc -l)


peakfile=$2

cut -f 1,2,3 ${peakfile}|coverageBed -split -b ${IP_sample} -a - | awk -v OFS='\t' -v SIZE=$IP_mappable_read_count '{print $1,$2,$3,($4*1000000/SIZE)*1000/($3-$2)}' > $outfile


import sys
import numpy as np

chroms = {}
bininfo = {}
with open(sys.argv[1]) as infile2:
    for line in infile2:
        line = line.strip()
        line_s = line.split('\t')
        chroms[line_s[-1]] = line_s[0]
        info1 = '-'.join(line_s[0:3])
        bininfo[info1] = int(line_s[-1])

interaction = {}
chrom_sum = {}

with open(sys.argv[2]) as infile2:
    for line in infile2:
        line = line.strip()
        line_s = line.split('\t')
        locusinfo = '_'.join(line_s[0:2])
        locusinfo2 = '_'.join(line_s[1::-1])
        interaction[locusinfo] = line_s[2]
        interaction[locusinfo2] = line_s[2]

        if not chroms[line_s[0]] in chrom_sum:
            chrom_sum[chroms[line_s[0]]] = 0
        if chroms[line_s[0]] == chroms[line_s[1]]:
            chrom_sum[chroms[line_s[0]]] += float(line_s[2])


sizeinfo = {}
with open("/nethome/yezhang_zhu/tools/ucsc_tools/mm9.chrom.sizes") as infile:
    for line in infile:
        line_s = line.strip().split("\t")
        sizeinfo[line_s[0]] = int(line_s[1])


TADs = []
with open(sys.argv[3]) as infile:
	for line in infile:
		line = line.strip().split('\t')
		TADs.append(line)

binsize = int(sys.argv[4])

for i in range(len(TADs)):
	for x in range(i+1,len(TADs),1):
	
		line_s = TADs[i]
		line2_s = TADs[x]
		if line_s[0] != line2_s[0]:
			continue

		if int(line_s[1]) % int(binsize) == 0:
			pass
		else:
			line_s[1] = str(int(binsize) * (int(line_s[1]) // int(binsize)))

		if int(line_s[2]) % int(binsize) == 0:
			pass
		else:
			line_s[2] = str(int(binsize) * (int(line_s[2]) // int(binsize)))

		if int(line_s[1])+binsize <= sizeinfo[line_s[0]]:
			chrom_start_info = line_s[0] + '-' + line_s[1] + '-' + str(int(line_s[1])+binsize)
		else:
			chrom_start_info = line_s[0] + '-' + line_s[1] + '-' + str(sizeinfo[line_s[0]])
		if int(line_s[2])+binsize <= sizeinfo[line_s[0]]:
			chrom_end_info = line_s[0] + '-' +  str(int(line_s[2])) + '-' + str(int(line_s[2])+binsize)
		else:
			chrom_end_info = line_s[0] + '-' +  str(int(line_s[2])) + '-' + str(sizeinfo[line_s[0]])

		bin_start = bininfo[chrom_start_info]
		bin_end = bininfo[chrom_end_info]

		if int(line2_s[1]) % int(binsize) == 0:
			pass
		else:
			line2_s[1] = str( int(int(binsize) * (int(line2_s[1]) // int(binsize))))

		if int(line2_s[2]) % int(binsize) == 0:
			pass
		else:
			line2_s[2] = str( int(int(binsize) * (int(line2_s[2]) // int(binsize))))
		
		if int(line2_s[1])+binsize <= sizeinfo[line2_s[0]]:
			chrom_start_info = line2_s[0] + '-' + line2_s[1] + '-' + str(int(line2_s[1])+binsize)
		else:
			chrom_start_info = line2_s[0] + '-' + line2_s[1] + '-' + str(sizeinfo[line2_s[0]])
		if int(line2_s[2])+binsize <= sizeinfo[line2_s[0]]:
			chrom_end_info = line2_s[0] + '-' +  str(int(line2_s[2])) + '-' + str(int(line2_s[2])+binsize)
		else:
			chrom_end_info = line2_s[0] + '-' +  str(int(line2_s[2])) + '-' + str(sizeinfo[line2_s[0]])

		bin_start2 = bininfo[chrom_start_info]
		bin_end2 = bininfo[chrom_end_info]

		sums = 0

		for m in range(bin_start,bin_end+1):
			for j in range(bin_start2,bin_end2+1):
				strs = str(m) + '_' + str(j)
				if strs in interaction:
					s1 = float(interaction[strs])
				else:
					s1 = 0
				sums += s1
		sums_norm = float(sums)/chrom_sum[line_s[0]]

		distance = abs((0.5 * (int(line_s[1]) + int(line_s[2]))) - (0.5 * (int(line2_s[1]) + int(line2_s[2]))))

		print("%s\t%s"%(sums_norm,distance))



import sys
import numpy as np
import math

chroms = {}
bininfo = {}
with open(sys.argv[1]) as infile2:
	for line in infile2:
		line = line.strip()
		line_s = line.split('\t')
		chroms[line_s[-1]] = line_s[0]
		info1 = '-'.join(line_s[0:3])
		bininfo[info1] = int(line_s[-1])

distance_x = {}
interaction = {}
with open(sys.argv[2]) as infile2:
    for line in infile2:
        line = line.strip()
        line_s = line.split('\t')
        locusinfo = '_'.join(line_s[0:2])
        locusinfo2 = '_'.join(line_s[1::-1])
        interaction[locusinfo] = line_s[2]
        interaction[locusinfo2] = line_s[2]

        froms = line_s[0]
        tos = line_s[1]
        if chroms[froms] == chroms[tos]:
            chrom = chroms[froms]
            distance = abs(int(froms) - int(tos))
            if not chrom in distance_x:
                distance_x[chrom] = {}

            if distance in distance_x[chrom]:
                distance_x[chrom][distance].append(float(line_s[2]))
            else:
                distance_x[chrom][distance] = [float(line_s[2])]

distance_x2 = {}
for chrom,value in distance_x.items():
    distance_x2[chrom] = {}
    for distance,value2 in value.items():
        distance_x2[chrom][distance] = np.mean(value2)


sizeinfo = {}
with open("/nethome/yezhang_zhu/tools/ucsc_tools/mm9.chrom.sizes") as infile:
	for line in infile:
		line_s = line.strip().split("\t")
		sizeinfo[line_s[0]] = int(line_s[1])

binsize = 20000

with open(sys.argv[3]) as infile:
	for line in infile:
		line_s = line.strip().split('\t')

		chrom = line_s[0]

		if int(line_s[1]) % int(binsize) == 0:
			pass
		else:
			line_s[1] = str(int(binsize) * (int(line_s[1]) // int(binsize)))

		if int(line_s[1])+binsize <= sizeinfo[line_s[0]]:
			chrom_start_info = line_s[0] + '-' + line_s[1] + '-' + str(int(line_s[1])+binsize)
		else:
			chrom_start_info = line_s[0] + '-' + line_s[1] + '-' + str(sizeinfo[line_s[0]])

		bins = bininfo[chrom_start_info]


		if int(line_s[4]) % int(binsize) == 0:
			pass
		else:
			line_s[4] = str(int(binsize) * (int(line_s[4]) // int(binsize)))

		if int(line_s[4])+binsize <= sizeinfo[line_s[0]]:
			chrom_start_info = line_s[0] + '-' + line_s[4] + '-' + str(int(line_s[4])+binsize)
		else:
			chrom_start_info = line_s[0] + '-' + line_s[4] + '-' + str(sizeinfo[line_s[0]])
		
		bins2 = bininfo[chrom_start_info]

		matrix2 = np.full((10 , 10 ), np.nan)

		for i in range(bins-5,bins+5):
			for j in range(bins2-5,bins2+5):
				distance = abs(i - j)
				index_i = i - bins
				index_j = j - bins2
				strs = str(i) + '_' + str(j)
				if strs in interaction:
					if distance in distance_x2[chrom]:
						matrix2[index_i,index_j] = float(interaction[strs]) /float(distance_x2[chrom][distance])
					else:
						matrix2[index_i,index_j] = float(interaction[strs])
		s1 = np.nanmean(matrix2)

		print(s1)



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


bed_lens = 0

with open(sys.argv[3]) as infile:
	bed_lens = len([i for i in infile.read().split('\n') if i != ""])

matrix = np.full((bed_lens, 50 , 50 ), np.nan)

binsize = 20000

with open(sys.argv[3]) as infile:
	index = 0
	for line in infile:
		line_s = line.strip().split('\t')
		chrom = line_s[0]


		if int(line_s[2]) % int(binsize) == 0:
			pass
		else:
			line_s[2] = str(int(binsize) * (int(line_s[2]) // int(binsize)))

		if int(line_s[2]) <= sizeinfo[line_s[0]]:
			chrom_end_info = line_s[0] + '-' +  str(int(line_s[2])-binsize) + '-' + line_s[2]
		else:
			chrom_end_info = line_s[0] + '-' +  str(int(line_s[2])-binsize) + '-' + str(sizeinfo[line_s[0]])

		bin_end = bininfo[chrom_end_info]

		bins_before = bin_end - 25
		bins_after = bin_end + 25

		bin_length = bins_after - bins_before 

		matrix2 = np.full((bin_length , bin_length ), np.nan)

		for i in range(bins_before,bins_after):
			for j in range(bins_before,bins_after):
				distance = abs(i - j)
				index_i = i - bins_before
				index_j = j - bins_before
				strs = str(i) + '_' + str(j)

				if strs in interaction:
					if distance in distance_x2[chrom]:
						matrix2[index_i,index_j] = float(interaction[strs]) /float(distance_x2[chrom][distance])
					else:
						matrix2[index_i,index_j] = float(interaction[strs])

		matrix[index,] = matrix2
		index += 1

matrix2 = np.nanmean(matrix,0)
np.savetxt(sys.argv[4],matrix2)


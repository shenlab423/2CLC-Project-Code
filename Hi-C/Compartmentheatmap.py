import sys
import numpy as np

def split_integer(m, n):
	assert n > 0
	quotient = int(m / n)
	remainder = m % n
	if remainder > 0:
		return [quotient] * (n - remainder) + [quotient + 1] * remainder
	if remainder < 0:
		return [quotient - 1] * -remainder + [quotient] * (n + remainder)
	return [quotient] * n

def meanArray(array):
	matrixout = np.full((50 , 50 ), np.nan)
	if len(array[0]) >= 50:
		splits = split_integer(len(array[0]),50)
		i = 0
		for m in range(50):
			j = 0
			for n in range(50):
				matrix2 = array[i:i+splits[m],j:j+splits[n]]
				s = np.nanmean(np.reshape(matrix2,(matrix2.size,)))
				matrixout[m,n] = s

				j += splits[n]
			i += splits[m]
		return(matrixout)
	else:
		return(matrixout)

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


matrix = np.full((22, 50 , 50 ), np.nan)


compartment_chr = {}
with open(sys.argv[3]) as infile:
	for line in infile:
		chrom = line.strip().split('\t')[0]
		if chrom == "chr1":
			if not chrom in compartment_chr:
				compartment_chr[chrom] = []
			compartment_chr[chrom].append('-'.join(line.strip().split('\t')[0:3]))

index = 0
for chrom,order in compartment_chr.items():
	bin_length = len(order)
	matrix2 = np.full((bin_length , bin_length ), np.nan)
	bins_order = [ bininfo[compartment_chr[chrom][i]] for i in range(len(compartment_chr[chrom])) ]

	
	index_i = 0

	for i in bins_order:
		index_j = 0
		for j in bins_order:
			distance = abs(i - j)
			strs = str(i) + '_' + str(j)
			if strs in interaction:
				if distance in distance_x2[chrom]:
					matrix2[index_i,index_j] = float(interaction[strs]) /float(distance_x2[chrom][distance])
				else:
					matrix2[index_i,index_j] = float(interaction[strs])
			index_j += 1
		index_i += 1


	matrix3 = meanArray(matrix2)
	matrix[index,] = matrix3
	index += 1

print(matrix)

matrix2 = np.nanmean(matrix,0)
np.savetxt(sys.argv[4],matrix2)


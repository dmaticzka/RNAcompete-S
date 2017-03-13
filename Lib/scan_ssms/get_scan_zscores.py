#!/home/hugheslab1/kate/bin/python

import sys
import csv

infile = sys.argv[1]
mean = float(sys.argv[2])
std = float(sys.argv[3])
cutoff = float(sys.argv[4])

with open(infile) as f:
	reader = csv.reader(f, delimiter='\t')
	for row in reader:
		binarized = []
		for base in row:
			if (float(base)-mean)/std > cutoff:
				binarized.append(1)
			else:
				binarized.append(0)
		print "\t".join(map(str,binarized))

f.close()

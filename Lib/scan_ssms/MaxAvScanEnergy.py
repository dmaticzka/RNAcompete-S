#!/usr/bin/python

import sys
from math import log

e = 2.718281828
PSEUDO = 0.00001
mu = 0 #This can be tuned 

#sys.argv[1] is the CisBp PFM
#sys.argv[2] is the Tab-Fasta 
#sys.argv[3] is the output name

#Read in PFM And Convert to Energy Matrix
def PFM2ematrix(filename):
	ematrix = []
	with open(filename, 'r') as pfm:
		for line in pfm:
			if line[0] != 'P':
				maxl = -1.0
				line = line.strip().split('\t')[1:]
				line = map(float, line)
				
				#find max frequency
				if max(line) > maxl:
					maxl = max(line) + PSEUDO
				
				#convert frequencies to energy scores
				A = -log(line[0]+ PSEUDO) + log(maxl)
				C = -log(line[1]+ PSEUDO) + log(maxl)
				G = -log(line[2]+ PSEUDO) + log(maxl)
				U = -log(line[3]+ PSEUDO) + log(maxl)
			
				#Add to matrix 
				ematrix.append([A,C,G,U])
	return(ematrix)

def MaxAvEnergyScanner(energymatrix, fastafilehandle, outputname, Motifprefix):
	with open(fastafilehandle,'r') as infasta:
		with open(outputname,'w') as Output:
			Header = ['Gene', (Motifprefix + '_MaxScore'),(Motifprefix + '_AverageScore'), (Motifprefix + '_SumScore')]
			Output.write('\t'.join(Header))
			Output.write('\n')
			BaseIndex = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
			for line in infasta:
				line = line.strip().split('\t')
				name = line[0]
				seq = line[1].upper()		
				pos = 0
				
				FwdScores = []
				
				#Score Forward Strand
				while pos < (len(seq) - len(energymatrix) + 1):
					Ei = 0
					subseq = seq[pos:(pos + len(energymatrix))]
					for subp in range(0, len(subseq)):
						base = subseq[subp]
						Ei += energymatrix[subp][BaseIndex[base]]
					subseqscore = 1/(1+(e**(Ei-mu)))
					FwdScores.append(subseqscore)
					pos += 1

				#Output Info
				
				outline = [name, seq]
				Output.write('\t'.join(map(str, outline)))
				Output.write('\t')
				Output.write(','.join(map(str, FwdScores)))
				Output.write('\n')


#Run the whole thing
MName = sys.argv[1].rstrip('.txt')
ematrix = PFM2ematrix(sys.argv[1])
MaxAvEnergyScanner(ematrix, sys.argv[2], sys.argv[3], MName)

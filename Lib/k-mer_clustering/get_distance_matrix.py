#!/home/hugheslab1/kate/bin/python

from numpy import *

import editdist
from util import *

import os,sys

# clusters k-mers
#
# fname: filename for list of k-mers
# nclusters: number of clusters
# use_kmer_scores: if True, then enrichment scores will be used for clustering
# clustAlgo: clustering method (to cluster eigenvectors), 'kmeans' or 'gmm' or 'dbscan'
#
# outputs an array of cluster labels for the k-mers
def run(fname,outfile_dist):
	kmers=loadtxt(fname,dtype='string')
	m = kmers.shape[0] # no. of kmers
	wts = array(kmers[:,2],float) # enrichment scores
	
	# compute pairwise distance matrix (where distance = the sum of the sequence and structure distance)
	D = zeros((m,m))
	for i in arange(m):
		print(i,m)
		for j in arange(i+1,m):
			#D[i,j]=editdist.distance(kmers[i,0],kmers[j,0])
			D[i,j]=compute_pfm_distance(kmers[i,0],kmers[j,0]) + compute_pfm_distance(kmers[i,1],kmers[j,1])
			D[j,i]=D[i,j]
	
	savetxt(outfile_dist,D)


if __name__ == "__main__":
	infile = sys.argv[1]
	outfile = sys.argv[2]
	run(infile,outfile)
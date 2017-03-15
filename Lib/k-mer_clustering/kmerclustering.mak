ID = $(shell head -2 ../info.tab | tail -1 | cut -f 2)
MODEL_FILE = model.tab
NUM_KMERS_TO_CLUSTER = 200 #$(shell cut -f 3 model.tab | grep "^0" | wc -l)
MAX_CLUSTERS = $(shell head -2 ../info.tab | tail -1 | cut -f 7)

targets  =  model.tab

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp) $(wildcard *.cache)

a:
	echo $(ID)
	echo $(MODEL_FILE)
	echo $(NUM_KMERS_TO_CLUSTER)

distance_matrix.txt: $(MODEL_FILE)
	cat $< \
	| head -$(NUM_KMERS_TO_CLUSTER) \
	| tr '\t' ' ' \
	> k.tmp;
	\
	python ../Lib/k-mer_clustering/get_distance_matrix.py k.tmp $@
	\
	rm -f k.tmp


plots: $(MODEL_FILE) distance_matrix.txt
	cat $< \
	| head -$(NUM_KMERS_TO_CLUSTER) \
	| tr '\t' ' ' \
	> k.tmp; \
	echo $(MAX_CLUSTERS); \
	../Lib/k-mer_clustering/clusterAndPlot.R k.tmp $(ID) $(MAX_CLUSTERS)

# clusterAndPlot.R produces labels.txt used here
clusters.txt: $(MODEL_FILE)
	cat $< \
	| head -$(NUM_KMERS_TO_CLUSTER) \
	| tr '\t' ' ' \
	> k.tmp;
	\
	paste k.tmp labels.txt \
	| tr ' ' '\t' \
	| sort -n -k 4,4 -k 3,3r \
	| ../Lib/perl_utilities/cap.pl kmer_seq,kmer_struct,feature_weight,cluster_ID \
	> $@;
	\
	mv clusters.txt $(ID)_clusters.txt;
	rm -f k.tmp labels.txt feat.txt

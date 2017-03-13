CLUST_ID = $(shell basename $(CURDIR))
ID = $(shell head -2 ../../info.tab | tail -1 | cut -f 2)

CLUST_FILE = ../$(ID)_clusters.txt

targets = seq_aligned.txt seq_aligned_trimmed.txt struct_aligned_trimmed.txt pfm_seq.txt pfm_struct.txt

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp) $(wildcard *.cache)

seq_aligned.txt: $(CLUST_FILE)
	cat $(CLUST_FILE) \
	| sed -e 1d \
	| select.pl -k 4 -eq $(CLUST_ID) \
	| cut -f 1 \
	| tr T U \
	> kmers_seq.tmp;
	\
	cat $(CLUST_FILE) \
	| sed -e 1d \
	| select.pl -k 4 -eq $(CLUST_ID) \
	| cut -f 2 \
	> kmers_struct.tmp;
	\
	python ../../Lib/cluster_ssms/align_seqstruct_kmers.py kmers_seq.tmp kmers_struct.tmp seq_aligned.txt struct_aligned.txt; \
	rm -f kmers_seq.tmp kmers_struct.tmp; \

seq_aligned_trimmed.txt: seq_aligned.txt
	python ../../Lib/cluster_ssms/trim_aligned.py $< $@;

struct_aligned_trimmed.txt: seq_aligned.txt
	python ../../Lib/cluster_ssms/trim_aligned.py struct_aligned.txt $@;

pfm_seq.txt: seq_aligned_trimmed.txt $(CLUST_FILE)
	cat $(CLUST_FILE) \
	| sed -e 1d \
	| select.pl -k 4 -eq $(CLUST_ID) \
	| cut -f 3 \
	> seq_kmers_weights.tmp;
	\
	cat seq_aligned_trimmed.txt \
	| paste - seq_kmers_weights.tmp \
	> seq_kmers_with_weights.tmp;
	\
	../../Lib/cluster_ssms/pfm_from_aligned_weighted.pl seq_kmers_with_weights.tmp \
	| lin.pl \
	| cap.pl PO,A,C,G,U \
	> $@;
	rm -f seq_aligned_trimmed.txt seq_kmers_weights.tmp seq_kmers_with_weights.tmp

pfm_struct.txt: struct_aligned_trimmed.txt $(CLUST_FILE)
	cat $(CLUST_FILE) \
	| sed -e 1d \
	| select.pl -k 4 -eq $(CLUST_ID) \
	| cut -f 3 \
	> struct_kmers_weights.tmp;
	\
	cat struct_aligned_trimmed.txt \
	| tr T U \
	| paste - struct_kmers_weights.tmp  \
	> struct_kmers_with_weights.tmp;
	\
	../../Lib/cluster_ssms/struct_pfm_from_aligned_weighted.pl struct_kmers_with_weights.tmp \
	| lin.pl \
	| cap.pl PO,B,E,H,L,M,R,T \
	> $@;
	rm -f seq_aligned.txt struct_aligned.txt struct_aligned_trimmed.txt struct_kmers_weights.tmp struct_kmers_with_weights.tmp

CLUST_ID = $(shell basename $(CURDIR))
ID = $(shell head -2 ../../info.tab | tail -1 | cut -f 2)

CLUST_IDS = $(shell cut -f 4 ../$(ID)_clusters.txt | sed 1d | uniq)

PFM_SEQ_FILE = pfm_seq.txt
PFM_STRUCT_FILE = pfm_struct.txt

VALID_SEQ_FILE = ../RBP_secondary_structure.txt.gz

NUM_SEQS = 10000
ZSCORE_CUTOFF = 2

targets =

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp)

a:
	echo $(CLUST_ID)
	echo $(ID)
	echo $(PFM_SEQ_FILE)
	echo $(PFM_STRUCT_FILE)
	echo $(VALID_SEQ_FILE)


scan_output_seq.tab: $(VALID_SEQ_FILE) $(PFM_SEQ_FILE)
	zcat $< \
	| cut -f 1,2 \
	| head -$(NUM_SEQS) \
	| tr T U \
	> input.tmp;
	python ../../Lib/scan_ssms/MaxAvScanEnergy.py $(PFM_SEQ_FILE) input.tmp $@;
	rm -f input.tmp;

scan_output_struct.tab: $(VALID_SEQ_FILE) $(PFM_STRUCT_FILE)
	zcat $< \
	| cut -f 1,3 \
	| head -$(NUM_SEQS) \
	> input.tmp;
	python ../../Lib/scan_ssms/MaxAvScanEnergy_struct.py $(PFM_STRUCT_FILE) input.tmp $@;
	rm -f input.tmp;

scan_output_summed.tab: scan_output_seq.tab scan_output_struct.tab
	cat scan_output_seq.tab \
	| sed -e 1d \
	| cut -f 3 \
	| tr , '\t' \
	> seq.tmp;
	\
	cat scan_output_struct.tab \
	| sed -e 1d \
	| cut -f 3 \
	| tr , '\t' \
	> struct.tmp;
	\
	python ../../Lib/scan_ssms/sum_matrix.py seq.tmp struct.tmp summed.tmp;
	\
	cat scan_output_seq.tab \
	| sed -e 1d \
	| cut -f 1 \
	| paste - summed.tmp\
	> $@;
	\
	rm -f seq.tmp struct.tmp summed.tmp;

max.txt: scan_output_summed.tab
	cat scan_output_summed.tab \
	| ../../Lib/perl_utilities/row_stats.pl -h 0 -argmax \
	| sed -e 1d \
	> $@;

diffs.txt: max.txt
	$(foreach c, $(CLUST_IDS), \
	   paste max.txt ../$(c)/max.txt \
	   | cut -f 2,4 \
	   | perl -ne 'chomp; @tabs=split; print $$tabs[1]-$$tabs[0]."\n";' \
	   | ../../Lib/perl_utilities/cap.pl diff_$(CLUST_ID)_$(c) \
	   > diff_$(c).tmp; \
	)
	paste diff_*.tmp > $@;
	rm -f diff_*.tmp scan_output_summed.tab scan_output_seq.tab scan_output_struct.tab;

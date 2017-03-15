ID = $(shell head -2 ../info.tab | tail -1 | cut -f 2)
CHILDREN = $(shell cut -f 4 $(ID)_clusters.txt | sed 1d | uniq)

targets = pfm_seq.txt pfm_struct.txt

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp) $(wildcard *.cache)


doit:
	$(foreach c, $(CHILDREN), \
	   cd $(c); \
	   echo $(c); \
	   make pfm_seq.txt; \
	   make pfm_struct.txt; \
	   cd ..; \
	)



maker:
	$(foreach c, $(CHILDREN), \
	   mkdir -p $(c); \
	   cd $(c); \
	   ln -sf ../../Lib/cluster_ssms/child.clust.mak Makefile; \
	   cd ..; \
	)

clust_sizes.txt:
	rm -f clust_sizes.tmp;
	$(foreach c, $(CHILDREN), \
	   wc -l $(c)/kmers_seq.tmp \
	   | cut -f 1 -d ' ' \
	   | ../Lib/perl_utilities/add_column.pl - -b -s '$(c)' \
	   >> clust_sizes.tmp; \
	)
	cat clust_sizes.tmp \
	| sort -nrk 2 \
	| ../Lib/perl_utilities/cap.pl clust_id,num_kmers \
	> $@;
	rm -f clust_sizes.tmp;
	


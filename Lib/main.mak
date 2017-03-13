ID = $(shell head -2 info.tab | tail -1 | cut -f 2)

annotated_secondary_structures:
	echo $(ID); \
	mkdir $(ID); \
	cd $(ID); \
	ln -sf ../Lib/structure_prediction_annotation/structure.mak Makefile; \
	make RBP_secondary_structure.txt.gz; \
	make POOL_secondary_structure.txt.gz; \
	cd ..; \

k-mer_model:
	echo $(ID); \
	cd $(ID); \
	rm -f Makefile; \
	ln -sf ../Lib/k-mer_model/kmermodel.mak Makefile; \
	make validation_auroc_seqstruct.txt; \
	cd ..; \

k-mer_clusters:
	echo $(ID); \
	cd $(ID); \
	rm -f Makefile; \
	ln -sf ../Lib/k-mer_clustering/kmerclustering.mak Makefile; \
	make plots; \
	make clusters.txt; \
	cd ..; \

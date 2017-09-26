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

cluster_ssms:
	echo $(ID); \
	cd $(ID); \
	rm -f Makefile; \
	ln -sf ../Lib/cluster_ssms/clusterssms.mak Makefile; \
	make maker; \
	make doit; \
	cd ..; \

scan_ssms:
	echo $(ID); \
	cd $(ID); \
	rm -f Makefile; \
	ln -sf ../Lib/scan_ssms/scanssms.mak Makefile; \
	make maker; \
	make doit; \
	make dist_histograms.pdf -B; \
	cd ..; \

final_ssm_assembly:
	echo $(ID); \
	cd $(ID); \
	rm -f Makefile; \
	ln -sf ../Lib/final_ssm_assembly/final_ssm_assembly.mak Makefile; \
	make pfm_assembled_seq_normalized.txt; \
	make pfm_assembled_struct_normalized.txt; \
	cd ..; \

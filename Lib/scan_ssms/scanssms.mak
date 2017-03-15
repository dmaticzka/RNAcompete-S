ID = $(shell head -2 ../info.tab | tail -1 | cut -f 2)
CHILDREN = $(shell cut -f 4 $(ID)_clusters.txt | sed 1d | uniq)

targets = all_diffs.txt dist_histograms.pdf

all: $(targets)

clean:
	rm -f $(targets) $(wildcard *.tmp) $(wildcard *.cache)

a:
	echo $(CHILDREN)

doit:
	$(foreach c, $(CHILDREN), \
	   cd $(c); \
	   echo $(c); \
	   make max.txt; \
	   cd ..; \
	) \
	$(foreach c, $(CHILDREN), \
	   cd $(c); \
	   echo $(c); \
	   make diffs.txt; \
	   cd ..; \
	) \


maker:
	$(foreach c, $(CHILDREN), \
	   cd $(c); \
	   rm -rf Makefile; \
	   ln -sf ../../Lib/scan_ssms/child.clust.mak Makefile; \
	   cd ..; \
	)

all_diffs.txt:
	paste */diffs.txt > $@;

dist_histograms.pdf: all_diffs.txt
	../Lib/scan_ssms/melt.R $< melted.txt;
	cat melted.txt \
	| sed -e 1d \
	| grep -v 1_1 \
	| grep -v 2_2 \
	| grep -v 3_3 \
	| grep -v 4_4 \
	| grep -v 5_5 \
	| sed 's/diff_//' \
	| sed 's/_/_other:/' \
	| tr _ '\t' \
	| ../Lib/perl_utilities/add_column.pl - -b -s 'key' \
	| sed 's/\t/:/' \
	| ../Lib/perl_utilities/cap.pl key_motif,other_motif,dist \
	> diff_input.txt;
	../Lib/scan_ssms/plot_dist_histograms.R diff_input.txt dist_histograms.png $@ all_dist_modes.txt;
	mv $@ $(ID)_dist_histograms.pdf; \
	rm -f melted.txt all_diffs.txt dist_histograms.png diff_input.txt;

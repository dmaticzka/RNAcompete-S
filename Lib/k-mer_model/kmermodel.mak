ID = $(shell head -2 ../info.tab | tail -1 | cut -f 2)
# define file with RBP selected sequences and file with background pool sequences
RBP_FILE = RBP_secondary_structure.txt.gz
POOL_FILE = POOL_secondary_structure.txt.gz

NUM_SEQS = $(shell zcat $(RBP_FILE) | wc -l)
NUM_VALID_SEQS = $(shell echo $$(( $(NUM_SEQS) / 3 )))
NUM_TRAIN_SEQS = $(shell echo $$(( $(NUM_SEQS) - $(NUM_VALID_SEQS) )))

REG_LAMBDA_L1 = $(shell head -2 ../info.tab | tail -1 | cut -f 5)
REG_LAMBDA_L2 = $(shell head -2 ../info.tab | tail -1 | cut -f 6)

shuffled_RBP_seqs.tab:
	zcat $(RBP_FILE) | shuf > $@

shuffled_pool_seqs.tab:
	zcat $(POOL_FILE) | shuf > $@

# train_RBP_vw_input.txt: shuffled_RBP_seqs.tab
# 	head -n $(NUM_TRAIN_SEQS) $< \
# 	| ../Lib/k-mer_model/make_seq_struct_vw_input.pl - 1 > $@

# train_pool_vw_input.txt: shuffled_pool_seqs.tab
# 	head -n $(NUM_TRAIN_SEQS) $< \
# 	| ../Lib/k-mer_model/make_seq_struct_vw_input.pl - -1 > $@

# train_vw_input.txt.gz: train_RBP_vw_input.txt train_pool_vw_input.txt
# 	cat $^ | gzip - > $@
# 	rm $^
# 	make valid_vw_input.txt.gz

train_vw_input.txt.gz: shuffled_RBP_seqs.tab shuffled_pool_seqs.tab
	( head -n $(NUM_TRAIN_SEQS) shuffled_RBP_seqs.tab \
	| ../Lib/k-mer_model/make_seq_struct_vw_input.pl - 1; \
	head -n $(NUM_TRAIN_SEQS) shuffled_pool_seqs.tab \
	| ../Lib/k-mer_model/make_seq_struct_vw_input.pl - -1 ) | gzip > $@
	make valid_vw_input.txt.gz

# valid_RBP_vw_input.txt: shuffled_RBP_seqs.tab
# 	tail -n $(NUM_VALID_SEQS) $< \
#         | ../Lib/k-mer_model/make_seq_struct_vw_input.pl - 1 > $@

# valid_pool_vw_input.txt: shuffled_pool_seqs.tab
# 	tail -n $(NUM_VALID_SEQS) $< \
#         | ../Lib/k-mer_model/make_seq_struct_vw_input.pl - -1 > $@

# valid_vw_input.txt.gz: valid_RBP_vw_input.txt valid_pool_vw_input.txt
# 	cat $^ | gzip - > $@
# 	rm $^
# 	rm shuffled_RBP_seqs.tab
# 	rm shuffled_pool_seqs.tab

valid_vw_input.txt.gz: shuffled_RBP_seqs.tab shuffled_pool_seqs.tab
	( tail -n $(NUM_VALID_SEQS) shuffled_RBP_seqs.tab \
        | ../Lib/k-mer_model/make_seq_struct_vw_input.pl - 1; \
	tail -n $(NUM_VALID_SEQS) shuffled_pool_seqs.tab \
	    | ../Lib/k-mer_model/make_seq_struct_vw_input.pl - -1 ) | gzip > $@
	rm shuffled_RBP_seqs.tab
	rm shuffled_pool_seqs.tab

predictor.l1l2.seqstruct: train_vw_input.txt.gz
	vw -d $< -b 30 -c  --passes 20 -f $@ --loss_function logistic  --compressed --keep q --keep r --keep t --keep s --keep u --keep v --keep w --l1 $(REG_LAMBDA_L1)  --l2 $(REG_LAMBDA_L2);

validation_predictions_seqstruct.txt: predictor.l1l2.seqstruct
	vw -d valid_vw_input.txt.gz -t -i $< -p $@ --compressed --invert_hash model_seqstruct.txt --keep q --keep r --keep t --keep s --keep u --keep v --keep w
	gzip model_seqstruct.txt;

validation_auroc_seqstruct.txt: validation_predictions_seqstruct.txt model.tab
	zcat valid_vw_input.txt.gz \
	| cut -f 1 -d ' ' \
	| paste - validation_predictions_seqstruct.txt \
	> data.tmp;
	../Lib/k-mer_model/get_AUROC.R data.tmp validation_auroc_seqstruct.txt $(ID)_seqstruct_ROC.pdf;
	rm -f data.tmp
	rm -f validation_predictions_seqstruct.txt

# to get a better formatted output model file
model.tab:
	zcat model_seqstruct.txt.gz \
	| tr : '\t' \
	| tr _ '\t' \
	| tr ^ '\t' \
	| tail -n +13 \
	| sort -nrk 5 \
	| cut -f 2,3,5 \
	> $@;
	rm -f model_seqstruct.txt.gz;
	rm -f *.cache
	rm -f predictor.l1l2.seqstruct

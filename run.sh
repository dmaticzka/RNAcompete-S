source /home/maticzkd/local_ncl/miniconda2/envs/rnacompete-s/bin/activate rnacompete-s

# RNAcompete-S reuquires libraries found in RNAcompete project
export PERL5LIB=~/projects/coop_rbns/RNAcompete/perl5/perl_tools:~/projects/coop_rnbs/RNAcompete/perl5/modules:~/scratch/rnacompete-s/RNAcompete/perl5/perl_tools:$PERL5LIB
# path to vw
export PATH=~/local/:$PATH

make -f Lib/main.mak
make -f Lib/main.mak k-mer_model
make -f Lib/main.mak k-mer_clusters
make -f Lib/main.mak cluster_ssms
make -f Lib/main.mak scan_ssms
make -f Lib/main.mak final_ssm_assembly

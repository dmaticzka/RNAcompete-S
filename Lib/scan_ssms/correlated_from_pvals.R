#!/usr/bin/Rscript

args = commandArgs(TRUE)

infile = args[1]
outfile = args[2]


library(Matrix)

df = read.table(infile, header=T)

motifs = levels(df$motif_ID_1)

combos = combn(motifs,2)

pvals = matrix(NA,length(motifs),length(motifs))

for (i in 1:ncol(combos) ) {
	df_subset = subset(df,motif_ID_1==combos[1,i]&motif_ID_2==combos[2,i])
	test = cor.test(df_subset$score_1,df_subset$score_2,alternative="greater")
	pval = test[["p.value"]]
	c1 = as.integer(unlist(strsplit(combos[1,i],'_',fixed=TRUE))[2])
	c2 = as.integer(unlist(strsplit(combos[2,i],'_',fixed=TRUE))[2])
	pvals[c1,c2] = pval
}

pvalsSym = forceSymmetric(pvals)

correlated = as.matrix(pvalsSym<0.01)

correlated2 = correlated
correlated2[correlated] <- 1

write.table(correlated2,outfile,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
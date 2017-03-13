#!/usr/bin/Rscript

args = commandArgs(TRUE)

library(gplots)

infile = args[1]
outfile = args[2]

df = read.table(infile)

png(outfile,width=400,height=800)

heatmap.2(as.matrix(df),
	Rowv=FALSE,
	Colv=FALSE,
	dendrogram="none",
	scale="none",
	trace="none",
	labRow=FALSE,
	labCol=FALSE,
	keysize=0.8,
	density.info="none")

dev.off()
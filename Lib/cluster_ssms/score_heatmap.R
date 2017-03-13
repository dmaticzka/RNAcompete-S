#!/usr/bin/Rscript

args = commandArgs(TRUE)

library(gplots)

infile = args[1]
outfile = args[2]

df = read.table(infile)

png(outfile,width=800,height=1000)

heatmap.2(as.matrix(df),Rowv=FALSE,Colv=FALSE,scale="none",trace="none",labRow=NULL)

dev.off()
#!/usr/bin/Rscript

args = commandArgs(TRUE)

df = read.table(args[1],header=T)

library(reshape2)

melted = melt(df)

write.table(melted,args[2],quote=F,sep="\t",row.names=F)
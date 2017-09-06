#!/usr/bin/env Rscript

args = commandArgs(TRUE)

df<-read.table(args[1])

outfile = args[3]

library(ROCR)

pred <- prediction(df$V2, df$V1)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")

pdf(outfile)
plot(perf)
abline(0,1,lty=2)
dev.off()

perf <- performance(pred, measure = "auc")

auc <- perf@y.values[[1]]

write(auc,file=args[2])

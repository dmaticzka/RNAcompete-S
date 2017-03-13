#!/usr/bin/Rscript

args = commandArgs(TRUE)

library(ggplot2)
library(plyr)

infile = args[1]
outfilepng = args[2]
outfilepdf = args[3]
outfilemodes = args[4]



df = read.table(infile,header=T)

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

means <- ddply(df, c("key_motif","other_motif"), summarise, position.mean=mean(dist,na.rm=T))
modes <- ddply(df, c("key_motif","other_motif"), summarise, position.mode=Mode(dist))
#medians <- ddply(df, c("key_motif","other_motif"), summarise, median=median(dist,na.rm=T))

write.table(modes,file=outfilemodes,quote=F,sep="\t",row.names=F)


pdf(outfilepdf)
g = ggplot(df,aes(x=dist)) +
	labs(x = "Distance between first base of key motif and first base other motif\n(key motif upstream of other motif)", y = "Count") +
	facet_grid(key_motif~other_motif)

g = g + geom_vline(aes(xintercept=0,linetype="solid"), size=1.5, colour="#CC0000")

g = g + geom_histogram()

plotdata = ggplot_build(g)$data[[1]]
modes$x = 15
modes$y = max(plotdata$count)-100

g = g + geom_vline(data=modes,
	aes(xintercept=position.mode,linetype="dotted"), size=1.5, colour="#3399FF" )

g = g + theme(legend.position="none")


g

dev.off()

png(outfilepng,width=800,height=800)

g = ggplot(df,aes(x=dist)) +
	geom_histogram() +
	labs(x = "Distance between first base of key motif and first base other motif\n(key motif upstream of other motif)", y = "Count") +
	facet_grid(key_motif~other_motif)
	
plotdata = ggplot_build(g)$data[[1]]

modes$x = 15
modes$y = max(plotdata$count)-100

g = g + geom_vline(data=modes,
	aes(xintercept=position.mode,linetype="dotted"), size=1.5, colour="#3399FF" )

g = g + geom_vline(aes(xintercept=0,linetype="solid"), size=1.5, colour="#CC0000")

g = g + theme(legend.position="none")

g = g + geom_text(data=modes,
		aes(label=paste("dist=",round(position.mode,2),sep=""),
			x=x,
			y=y),
		fontface="bold",
		inherit.aes=FALSE)

g

dev.off()
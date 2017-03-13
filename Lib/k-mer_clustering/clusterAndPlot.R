#!/usr/bin/Rscript

### DEPENDENCIES ####
# smacof, mclust and calibrate
#####################


args <- commandArgs(TRUE)

infname = args[1]
outfname = args[2]
maxclusts = args[3]
print(maxclusts)

# load all packages
library(mclust)
#library(calibrate)
# load smacof files (for weighted MDS)
files<-list.files('../Lib/k-mer_clustering/smacof')
for(f in files){source(paste('../Lib/k-mer_clustering/smacof/',f,sep=''))}


kmers<-read.table(infname) # load list of kmers and their enrichment scores
e<-kmers$V3/max(kmers$V3) #enrichment scores
dist<-read.table("distance_matrix.txt") # load (PFM) distance matrix

# compute embeddings using smacof method
res <- smacofSym(dist,ndim=10,weightmat=outer(e,e),verbose=FALSE,itmax = 5000)
feat<-res$conf

write.table(feat,"feat.txt",quote=F,sep="\t",row.names=F)

#call the weighted clustering (GMM) method
r<-maxclusts
model<-Mclust(feat,G=1:r,modelNames=c('EEI', 'VEI', 'EVI', 'VVI'),wts=e)
summary(model)

write(map(model$z),ncolumns=1,file=paste("labels.txt",sep="")) # write the labels to file

# add kmer text to plot
kmers$labs = paste(kmers$V1,kmers$V2,sep="_")

pdf(paste(outfname,'_mclust2Dplot','.pdf',sep=""))
mclust2Dplot(data=feat[,c(1,2)], parameters = model$parameters, classification = model$classification )

zz<-map(model$z)
cids<-unique(zz)
tids<-list()
for(id in cids){
	tids<-which(zz==id)[1:5]
	text(feat[tids,1],feat[tids,2], labels=array(kmers$labs[tids]), cex = 1)
}

title(main=paste(outfname," clusters, plotted on first 2 dimensions",sep=""),
  	xlab="Coordinate 1", ylab="Coordinate 2")

#mclust2Dplot(data=feat[,c(1,3)], parameters = model$parameters, classification = model$classification )
#mclust2Dplot(data=feat[,c(1,4)], parameters = model$parameters, classification = model$classification )
#mclust2Dplot(data=feat[,c(1,5)], parameters = model$parameters, classification = model$classification )
#mclust2Dplot(data=feat[,c(2,3)], parameters = model$parameters, classification = model$classification )
#mclust2Dplot(data=feat[,c(2,4)], parameters = model$parameters, classification = model$classification )
#mclust2Dplot(data=feat[,c(2,5)], parameters = model$parameters, classification = model$classification )
#mclust2Dplot(data=feat[,c(3,4)], parameters = model$parameters, classification = model$classification )
#mclust2Dplot(data=feat[,c(3,5)], parameters = model$parameters, classification = model$classification )
#mclust2Dplot(data=feat[,c(4,5)], parameters = model$parameters, classification = model$classification )
dev.off()

if(FALSE){
# 2D projecion for plotting
DR = MclustDR(model)

pdf(paste(outfname,'_DRplot','.pdf',sep=""))
plot(DR)
plot(DR, what = "evalues")
plot(DR, what = "pairs")
plot(DR, what = "contour")
plot(DR, what = "classification", ngrid = 200)
plot(DR, what = "boundaries", ngrid = 200)
plot(DR, what = "density")
plot(DR, what = "density", dimens = 2)
dev.off()



#pdf(paste(outfname,'_mclustplot','.pdf',sep=""))
#plot(model)
#dev.off()

# combi object to combine mixture components
combi <- clustCombi(feat, model)


object <- DR
x <- object$x
p <- ncol(x)
n <- nrow(x)
G <- object$G
y <- object$mixcomp


print('G=')
print(head(G))
print('p=')
print(head(p))
print('n=')
print(head(n))


class <- as.numeric(object$class)
nclass <- length(table(class))

dir <-object$dir
numdir <- object$numdir
dimens <- seq(numdir)

print(head(dir))
print(length	(dir))
print('number of clusters:')
print(nclass)
print('cluster labels:')
print(head(class))

dimens <- dimens[1:2]
dir <- object$dir[,dimens,drop=FALSE]
print(head(dir))
par <- projpar.MclustDR(object, dimens)
Mu <- par$mean
Sigma <- par$variance

# draw contours for each class or cluster
asp = NULL
symbols <- mclust.options('classPlotSymbols')
colors <- mclust.options('classPlotColors')



pdf(paste(outfname,'.pdf',sep=""))

plot.window(xlim = 1.2*range(dir[,1]), ylim = 1.2*range(dir[,2]), asp = asp)
plot(dir, type = "n", xlim = par("usr")[1:2], ylim = par("usr")[3:4], xaxs = "i", yaxs = "i", main=outfname)
#if(length(dimens)>1) # no contours for 1D projection...
#{
	for(k in seq(nclass))
	{
		i <- which(object$class2mixcomp == k)
		parameters <- list(pro = object$pro[i]/sum(object$pro[i]), mean = t(par$mean[i,,drop=FALSE]), variance = list(G = length(i), d = 2, sigma = par$variance[,,i,drop=FALSE]))
		surfacePlot(dir, parameters, col = gray(.7), nlevels = 5, xlim = par("usr")[1:2], ylim = par("usr")[3:4], asp = asp, add = TRUE)
	}
#}
points(dir, col = colors[class], pch = symbols[class],cex=2*e)

# add kmer text to plot
kmers$labs = paste(kmers$V1,kmers$V2,sep="_")

zz<-map(model$z)
cids<-unique(zz)
tids<-list()
rot<-0 # rotate text for 1D projection
if(length(dimens)<2){rot<-90}
for(id in cids){
	tids<-which(zz==id)[1:5]
	textxy(dir[tids,1],dir[tids,2], srt=rot, labs=array(kmers$labs[tids]), cex = 0.3,  m = c(0, 0))
}

dev.off()

}


#### script ends here

if(FALSE){
###### ENTROPY PLOTS AND STUFF ##################

# further combine/visualize the clusters
source('util.R')

z <- combi$MclustOutput$z
combiM <- combi$combiM
ent <- numeric()
Kmax <- ncol(z)
Kmax <- ncol(z)


abc <- "standard"
reg <-2 # default
if (Kmax == 2) reg = NULL # ??

# compute entropy
ent <- numeric()
z0 <- z
for (K in Kmax:1)
{
	z0 <- t(combiM[[K]] %*% t(z0))
	ent[K] <- -sum(xlog(z0))
}
# compute normalized entropy
if (any(abc == "normalized"))
{
	mergedn <-numeric()
	z0 <- z
	for (K in (Kmax-1):1)
	{
		z0 <- t(combiM[[K+1]] %*% t(z0))
		mergedn[K] = sum(sapply(map(z0),function(x) any(which(as.logical(combiM[[K]][rowSums(combiM[[K]])==2,]))==x)))
	}

}

if (any(abc == "standard"))
{
	plot(1:Kmax, ent, xlab = "Number of clusters", ylab = "Entropy")
	if (any(reg == 2))
	{
		pcwsreg <- pcws2_reg(1:Kmax,ent)
		lines(1:pcwsreg$c, pcwsreg$a1*(1:pcwsreg$c) + pcwsreg$b1, lty = 2, col = "red")
		lines(pcwsreg$c:Kmax, pcwsreg$a2*(pcwsreg$c:Kmax) + pcwsreg$b2, lty = 2, col = "red")
	}
	if (any(reg == 3))
	{
		pcwsreg <- pcws3_reg(1:Kmax,ent)
		lines(1:pcwsreg$c1, pcwsreg$a1*(1:pcwsreg$c1) + pcwsreg$b1, lty = 2, col = "blue")
		lines(pcwsreg$c1:pcwsreg$c2, pcwsreg$a2*(pcwsreg$c1:pcwsreg$c2) + pcwsreg$b2, lty = 2, col = "blue")
		lines(pcwsreg$c2:Kmax, pcwsreg$a3*(pcwsreg$c2:Kmax) + pcwsreg$b3, lty = 2, col = "blue")
	}
	#plot(1:(Kmax-1), ent[2:Kmax]-ent[1:(Kmax-1)], xlab = "Number of clusters", ylab = "Difference in entropy")
	title("Entropy plot", outer=TRUE)
}

if (any(abc == "normalized"))
{
	plot(cumsum(c(0,mergedn)), ent, xlab = "Cumul. count of merged obs.", ylab = "Entropy")
	if (any(reg == 2))
	{
		X <- cumsum(c(0,mergedn))
		pcwsreg <- pcws2_reg(X,ent)
		lines(X[1:pcwsreg$c], pcwsreg$a1*(X[1:pcwsreg$c]) + pcwsreg$b1, lty = 2, col = "red")
		lines(X[pcwsreg$c:Kmax], pcwsreg$a2*(X[pcwsreg$c:Kmax]) + pcwsreg$b2, lty = 2, col = "red")
	}
	if (any(reg == 3))
	{
		X <- cumsum(c(0,mergedn))
		pcwsreg <- pcws3_reg(X,ent)
		lines(X[1:pcwsreg$c1], pcwsreg$a1*(X[1:pcwsreg$c1]) + pcwsreg$b1, lty = 2, col = "blue")
		lines(X[pcwsreg$c1:pcwsreg$c2], pcwsreg$a2*(X[pcwsreg$c1:pcwsreg$c2]) + pcwsreg$b2, lty = 2, col = "blue")
		lines(X[pcwsreg$c2:Kmax], pcwsreg$a3*(X[pcwsreg$c2:Kmax]) + pcwsreg$b3, lty = 2, col = "blue")
	}
	#plot(1:(Kmax-1), (ent[2:Kmax]-ent[1:(Kmax-1)])/mergedn, xlab = "Number of clusters", ylab = "Normalized difference in entropy")
	title("Normalized entropy plot", outer=TRUE)
}

}

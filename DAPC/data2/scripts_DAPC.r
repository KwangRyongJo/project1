library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")

setwd("path")
datafile <- "CC_QC_SK_US_JP_AABB_imputed_393g_4440m.txt" 
dat <- read.table(datafile)
head(dat)

obj <- df2genind(dat, ploidy=4, sep="")
obj

X <- tab(obj, NA.method="mean")

# use find.clusters to identify clusters
grp <- find.clusters(obj, max.n.clust=40) # choose 200, then 6
##Choose the number PCs to retain (>= 1): 200
##Choose the number of clusters (>=2): 6


dapc1 <- dapc(X, grp$grp) # choose 120, then 6
##Choose the number PCs to retain (>=1): 120
##Choose the number discriminant functions to retain (>=1): 6


dapc1
scatter(dapc1)

# saving	
pdf("SK_393g_4440m.pdf", width=30, height=15) # Open a PDF
p <- scatter(dapc1, label.inds = list(air = 1))
print(p)
dev.off()

# SNP contributions
loadingplot(dapc1$var.contr)
loadingplot(tail(dapc1$var.contr, 100), main="Loading plot - last 100 SNPs")
write.csv(dapc1$var.contr, file = "loadings_393g_4440m.csv")

myCol <- c("blue","red","orange", "cyan","green","brown")
scatter(dapc1, label.inds = list(air = 1), posi.da="bottomright", bg="white",
        pch=17:22, col=myCol, scree.pca=TRUE,
        posi.pca="topright")


scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
        cex=3,clab=0, leg=TRUE, posi.leg="topright", txt.leg=paste("Cluster",1:3))	


scatter(dapc1, label.inds = list(air = 1), ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=myCol, solid=.4, cex=3, clab=0,
        mstree=TRUE, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg="bottomright", txt.leg=paste("Cluster",1:6))
par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=3, lwd=8, col="white")
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=myCol)

myInset <- function(){
    temp <- dapc1$pca.eig
    temp <- 100* cumsum(temp)/sum(temp)
    plot(temp, col=rep(c("black","lightgrey"),
                       c(dapc1$n.pca,1000)), ylim=c(0,100),
         xlab="PCA axis", ylab="Cumulated variance (%)",
         cex=1, pch=20, type="h", lwd=2)
}
add.scatter(myInset(), posi="topright",
            inset=c(-0.03,-0.01), ratio=.18,
            bg=transp("white"))


scatter(dapc1,1,1, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)


clusters <- round(dapc1$posterior,6)
write.csv(clusters, file="clusters.csv")

# which are the most 'admixed' individuals
# Let us consider as admixed individuals 
# having no more than 90% of probability of membership 
# in a single cluster:
temp <- which(apply(dapc1$posterior,1, function(e) all(e<0.9)))
temp

assignplot(dapc1, subset=1:50)

compoplot(dapc1, posi="bottomright",
          txt.leg=paste("Cluster", 1:6), lab="",
          ncol=1, xlab="potato clones")


compoplot(dapc1, subset=temp, posi="bottomright",
          txt.leg=paste("Cluster", 1:6),
          ncol=2)

assignplot(dapc1)


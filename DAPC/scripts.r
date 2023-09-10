library("ape")
library("pegas")
library("seqinr")
library("ggplot2")

library("adegenet")

setwd("D:/project1/DAPC")
datafile <- "QC_AATT_127g_6671m.txt" # cc2
dat <- read.table(datafile)
head(dat)

obj <- df2genind(dat, ploidy=4, sep="")
obj

# genind2df(obj, sep="|")

X <- tab(obj, NA.method="mean")

# use find.clusters to identify clusters
grp <- find.clusters(obj, max.n.clust=40) # choose 200, then 3
##Choose the number PCs to retain (>= 1): 200
##Choose the number of clusters (>=2): 3
##fig-1

dapc1 <- dapc(X, grp$grp) # choose 60, then 3
##Choose the number PCs to retain (>=1): 60
##Choose the number discriminant functions to retain (>=1): 3
##fig-2

dapc1
scatter(dapc1,label.inds = list(air = 1))

# saving	
pdf("DAPC_K.pdf", width=30, height=15) # Open a PDF
p <- scatter(dapc1, label.inds = list(air = 1))

print(p)
dev.off()

# SNP contributions
loadingplot(dapc1$var.contr)
loadingplot(tail(dapc1$var.contr, 100), main="Loading plot - last 100 SNPs")
#write.csv(dapc1$var.contr, file = "10-DAPC/cc.csv")

myCol <- c("blue","red","orange")
scatter(dapc1, label.inds = list(air = 1), posi.da="bottomleft", bg="white",
        pch=17:22, col=myCol, scree.pca=TRUE,
        posi.pca="bottomright")

scatter(dapc1, label.inds = list(air = 1), scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,
cex=3,clab=0, leg=TRUE, posi.leg="topleft", txt.leg=paste("Cluster",1:3))


scatter(dapc1, label.inds = list(air = 1), ratio.pca=0.3, bg="white", pch=20, cell=0,
        cstar=0, col=myCol, solid=.4, cex=3, clab=0,
        mstree=TRUE, scree.da=FALSE, posi.pca="bottomright",
        leg=TRUE, posi.leg="topleft", txt.leg=paste("Cluster",1:3))
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
add.scatter(myInset(), posi="bottomleft",
            inset=c(-0.03,-0.01), ratio=.18,
            bg=transp("white"))

scatter(dapc1,1,1, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)

summary(dapc1)

# which are the most 'admixed' individuals
# Let us consider as admixed individuals 
# having no more than 90% of probability of membership 
# in a single cluster:
temp <- which(apply(dapc1$posterior,1, function(e) all(e<0.9)))
temp

################################# 
library(adegenet)
library(poppr)
library(dplyr)
library(hierfstat)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(scales)

setwd("D:/project1/DAPC")
dat = read.csv("Nieuwe map/KAJ.csv")
str(dat)

# Create vectors of individual and site labels.
ind = as.character(dat$ID) # individual ID
site = as.character(dat$Site) # site ID
 
# Convert data.frame to genind object. 
kaj = df2genind(dat, ploidy = 4, ind.names = ind, pop = site, sep = "")
kaj$tab[1:5, 1:10]

x <- tab(kaj, NA.method="mean")
 
# Perform PCA
pca1 = dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)
 
# Analyse how much percent of genetic variance is explained by each axis
percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1))
 
 
# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca1$li)
 
# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")
 
# Add a column containing individuals
ind_coords$Ind = indNames(kaj)
 
# Add a column with the site IDs
ind_coords$Site = kaj$pop
# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)
 
 
# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
 
# Define colour palette
cols = brewer.pal(nPop(kaj), "Set1")

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")
 
# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15) 
)
 
# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
    geom_hline(yintercept = 0)+
    geom_vline(xintercept = 0)+
    # spider segments
    geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
    # points
    geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
    # centroids
    geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
    # colouring
    scale_fill_manual(values = cols)+
    scale_colour_manual(values = cols)+
    # custom labels
    labs(x = xlab, y = ylab)+
    ggtitle("Potato PCA")+
    # custom theme
    ggtheme


# Perform cross validation to find the optimal number of PCs to retain in DAPC
set.seed(123)
x = tab(kaj, NA.method = "mean")
crossval = xvalDapc(x, kaj$pop, result = "groupMean", xval.plot = TRUE)

# Number of PCs with best stats (lower score = better)
crossval$`Root Mean Squared Error by Number of PCs of PCA`
#        20        40        60        80       100       120       140       160       180       200       220 
# 0.2007724 0.1633198 0.1748439 0.1710502 0.1386689 0.1824984 0.1902544 0.1895499 0.2425877 0.3260219 0.4613567 

crossval$`Number of PCs Achieving Highest Mean Success`
# [1] "100"

crossval$`Number of PCs Achieving Lowest MSE`
# [1] "100"
numPCs = as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)

# Run a DAPC using site IDs as priors
dapc1 = dapc(kaj, kaj$pop, n.pca = numPCs, n.da = 3)

# Analyse how much percent of genetic variance is explained by each axis
percent = dapc1$eig/sum(dapc1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,60),
        names.arg = round(percent, 1))

# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(dapc1$ind.coord)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2")

# Add a column containing individuals
ind_coords$Ind = indNames(kaj)

# Add a column with the site IDs
ind_coords$Site = kaj$pop

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2) ~ Site, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))

# Define colour palette
cols = brewer.pal(nPop(kaj), "Set2")

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Site), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Site), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Site), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("DAPC analysis of potato cultivars bred in Korea, Japan and US")+
  # custom theme
  ggtheme









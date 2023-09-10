library(adegenet)
library(poppr)
library(dplyr) 
library(hierfstat) 
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(factoextra)

setwd("D:/project1/KLFDAPC/")
dat = read.csv("KAJE.csv")
# Create vectors of individual and site labels.
ind = as.character(dat$ID) # individual ID
site = as.character(dat$Site) # site ID
 
# Convert data.frame to genind object. 
kaje = df2genind(dat, ploidy = 4, ind.names = ind, pop = site, sep = "")
x <- tab(kaje, NA.method="mean")

# Perform PCA
pca1 = dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 10)
 
# Analyse how much percent of genetic variance is explained by each axis
percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1))
 
# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca1$li)

write.csv(ind_coords, file = "KAJE_PC1_10.csv")

library(KLFDAPC)
library(PCAviz)
library(ggplot2)
library(cowplot)

setwd("D:/project1/KLFDAPC/")
regmappcs1 <- read.csv(file = "KAJE_PC1_10x.csv")
str(regmappcs1)
##'data.frame':   272 obs. of  22 variables:
## $ X               : int  1 2 3 4 5 6 7 8 9 10 ...
## $ array_id        : int  641 642 643 644 645 646 647 648 649 650 ...
## $ ecotype_id      : chr  "AlpineRusset" "BeaconChipper" "BlazerRusset" "CentennialRusset" ...
## $ median_intensity: int  476 556 542 478 560 363 335 482 360 353 ...
## $ latitude        : num  55.8 55.7 63 63 57.7 ...
## $ longitude       : num  14.2 13.4 18.3 18.3 15 ...
## $ nativename      : chr  "Fri 1" "FlyA 3" "FaU 4" "FaL 1" ...
## $ firstname       : chr  "Alison" "Alison" "Alison" "Alison" ...
## $ surname         : chr  "Anastasio" "Anastasio" "Anastasio" "Anastasio" ...
## $ site            : chr  "Frisaboda" "Flyinge" "Faberget-Upper" "Faberget-Lower" ...
## $ region          : chr  "US" "US" "US" "US" ...
## $ country         : chr  "US" "US" "US" "US" ...
## $ PC1             : num  -9.36 5.56 -9.53 -5.19 -6.21 ...
## $ PC2             : num  -0.747 -3.127 -2.651 -1.999 -0.746 ...
## $ PC3             : num  13.18 6.79 15.33 7.39 5.8 ...
## $ PC4             : num  2.04 5.49 -1.02 -2.54 -1.56 ...
## $ PC5             : num  -6.885 2.446 -0.304 3.682 -4.722 ...
## $ PC6             : num  -1.938 -6.338 -3.68 -3.481 -0.882 ...
## $ PC7             : num  -12.13 2.54 -4.35 3.76 -9.93 ...
## $ PC8             : num  1.67 -2.76 8.11 3.2 -2.77 ...
## $ PC9             : num  -3.084 1.571 -5.589 0.193 8.611 ...
## $ PC10            : num  1.63 -11.26 3.57 6.19 -1.08 ...

regmappcs1$region <- as.factor(regmappcs1$region) # Convert character column to factor
regmappcs1$country <- as.factor(regmappcs1$country) # Convert character column to factor
regmappcs1$array_id <- as.character(regmappcs1$array_id) # Convert int column to character
regmappcs1$ecotype_id <- as.character(regmappcs1$ecotype_id) # Convert int column to character
str(regmappcs1) 

#head(regmappcs1)

### we have to remove the not well represented individuals, which one country only has one individuals.
table(regmappcs1$country)

regmapviz=pcaviz(dat = regmappcs1)

summary(regmapviz)

### 3 countries
regmapviz$data$country=factor(regmapviz$data$country,levels = unique(regmapviz$data$country))

regmapviz1=pcaviz(dat = regmapviz$data)

### The first 10 PCs, we will produce 3 genetic features for visualization
normalize <- function(x) {
return ((x - min(x)) / (max(x) - min(x)))
}


pcanorm=apply(regmapviz1$data[,13:22], 2, normalize)

reg_kmat <- kmatrixGauss(pcanorm,sigma=2)

reg_klfdapc=KLFDA(reg_kmat, y=regmapviz1$data$country, r=3, knn = 5)

### Visualization of the population structure
# project the first two genetic features onto the plane
regmapviz1$data$PC1=reg_klfdapc$Z[,1]

regmapviz1$data$PC2=reg_klfdapc$Z[,2]

plot(regmapviz1,coord=c("PC1","PC2"),group =NULL,draw.points =FALSE,label = "ecotype_id",plot.title=" KLFDAPC 1 vs. KLFDAPC 2")+xlab("KLFDAPC 1")+ylab("KLFDAPC 2")



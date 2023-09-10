setwd("path")
library(Hapi)
datafile1 <- "SNP_pos.txt"
data_c1 <- read.table(datafile1,header=T)
head(data_c1)

datafile2 <- "potato_chr.txt"
data_c2 <- read.table(datafile2,header=T)
head(data_c2)

### view gamete cells 
hapiGameteView(chr = data_c2, hap = data_c1)

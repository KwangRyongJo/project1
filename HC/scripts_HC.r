library(StAMPP)	
potato.freq <- stamppConvert("QC_CC_AABB_131g_4946m_imputed4.csv")
 
# calculate pair-wise Nei’s genetic distances between individuals (FALSE)
potato.ind.d <- stamppNeisD(potato.freq, pop=FALSE)
 
library("factoextra")
obj <- as.dist(potato.ind.d)


library(ape)
library(phangorn)
library(seqinr)
library(dendextend)
library(factoextra)
library(circlize)


dend <- obj %>% 
    hclust(method = "ward.D2") %>% as.dendrogram %>%
    set("branches_k_color", value = c("cyan", "blue", "orange"), k=3) %>% 
    set("branches_lwd", 2.2) %>%     
    set("labels_colors", c('brown', 'red', 'brown', 'brown', 'green', 'red', 'brown', 'red', 'brown', 'brown', 'red', 'red', 'brown', 'black', 'brown', 'brown', 'red', 'red', 'red', 'brown', 'red', 'red', 'red', 'red', 'cyan', 'black', 'black', 'brown', 'brown', 'brown', 'red', 'red', 'brown', 'red', 'green', 'black', 'green', 'red', 'purple', 'purple', 'purple', 'purple', 'black', 'black', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'black', 'brown', 'red', 'brown', 'brown', 'brown', 'black', 'red', 'brown', 'red', 'red', 'brown', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'black', 'purple', 'purple', 'black', 'black', 'black', 'brown', 'green', 'black', 'brown', 'black', 'black', 'cyan', 'green', 'black', 'brown', 'black', 'black', 'black', 'black', 'green', 'brown', 'black', 'red', 'brown', 'green', 'black', 'brown', 'red', 'brown', 'red', 'brown', 'black', 'brown', 'black', 'black', 'purple', 'brown', 'brown', 'brown', 'brown', 'green', 'brown', 'brown', 'green', 'brown', 'black', 'brown', 'brown', 'brown', 'brown'
)) %>% 
    set("labels_cex", c(0.75)) %>% 
    set("leaves_cex", c(1.25)) %>%
    set("leaves_pch", 19) %>% set("leaves_col", c('brown', 'red', 'brown', 'brown', 'green', 'red', 'brown', 'red', 'brown', 'brown', 'red', 'red', 'brown', 'black', 'brown', 'brown', 'red', 'red', 'red', 'brown', 'red', 'red', 'red', 'red', 'cyan', 'black', 'black', 'brown', 'brown', 'brown', 'red', 'red', 'brown', 'red', 'green', 'black', 'green', 'red', 'purple', 'purple', 'purple', 'purple', 'black', 'black', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'black', 'brown', 'red', 'brown', 'brown', 'brown', 'black', 'red', 'brown', 'red', 'red', 'brown', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'black', 'purple', 'purple', 'black', 'black', 'black', 'brown', 'green', 'black', 'brown', 'black', 'black', 'cyan', 'green', 'black', 'brown', 'black', 'black', 'black', 'black', 'green', 'brown', 'black', 'red', 'brown', 'green', 'black', 'brown', 'red', 'brown', 'red', 'brown', 'black', 'brown', 'black', 'black', 'purple', 'brown', 'brown', 'brown', 'brown', 'green', 'brown', 'brown', 'green', 'brown', 'black', 'brown', 'brown', 'brown', 'brown'
))


par(mar = rep(0, 4))
circlize_dendrogram(dend, dend_track_height = 0.8)

######################################## 
library(StAMPP)	
potato.freq <- stamppConvert("J_K_US_410g_4043m2.csv")
 
# calculate pair-wise Nei’s genetic distances between individuals (FALSE)
potato.ind.d <- stamppNeisD(potato.freq, pop=FALSE)
 
library("factoextra")
obj <- as.dist(potato.ind.d)


library(ape)
library(phangorn)
library(seqinr)
library(dendextend)
library(factoextra)
library(circlize)


dend <- obj %>% 
    hclust(method = "ward.D2") %>% as.dendrogram %>%
    set("branches_k_color", value = c("darkgrey", "blue", "orange", "brown", "green", "purple", "yellow", "light blue", "red", "pink"), k=10) %>% 
    set("branches_lwd", 2.0) %>%     
    set("labels_colors", c('red', 'black', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'black', 'brown', 'red', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'black', 'brown', 'brown', 'brown', 'black', 'brown', 'brown', 'red', 'green', 'red', 'brown', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'red', 'red', 'red', 'brown', 'red', 'red', 'brown', 'red', 'red', 'red', 'red', 'green', 'brown', 'brown', 'red', 'red', 'red', 'red', 'red', 'green', 'red', 'brown', 'red', 'red', 'black', 'red', 'black', 'red', 'black', 'red', 'brown', 'brown', 'brown', 'brown', 'red', 'red', 'brown', 'red', 'black', 'red', 'brown', 'brown', 'brown', 'cyan', 'red', 'brown', 'red', 'red', 'red', 'brown', 'red', 'red', 'purple', 'brown', 'red', 'red', 'red', 'red', 'brown', 'brown', 'brown', 'red', 'orange', 'orange', 'orange', 'orange', 'green', 'red', 'brown', 'red', 'green', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'brown', 'brown', 'brown', 'red', 'red', 'red', 'black', 'brown', 'orange', 'brown', 'brown', 'brown', 'brown', 'red', 'red', 'black', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'orange', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'brown', 'red', 'red', 'green', 'red', 'red', 'red', 'green', 'green', 'green', 'green', 'green', 'green', 'green', 'green', 'brown', 'orange', 'green', 'green', 'brown', 'black', 'brown', 'brown', 'green', 'green', 'green', 'green', 'green', 'green', 'green', 'brown', 'green', 'green', 'brown', 'green', 'green', 'green', 'green', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'red', 'orange', 'red', 'brown', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'black', 'black', 'brown', 'brown', 'brown', 'brown', 'brown', 'black', 'brown', 'brown', 'purple', 'black', 'red', 'black', 'green', 'brown', 'brown', 'brown', 'brown', 'cyan', 'brown', 'black', 'cyan', 'orange', 'brown', 'brown', 'green', 'green', 'red', 'red', 'black', 'brown', 'brown', 'red', 'brown', 'green', 'green', 'green', 'brown', 'brown', 'red', 'purple', 'red', 'brown', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'orange', 'purple', 'brown', 'brown', 'purple', 'orange', 'brown', 'cyan', 'brown', 'green', 'brown', 'brown', 'purple', 'purple', 'brown', 'purple', 'brown', 'purple', 'purple', 'purple', 'brown', 'black', 'brown', 'brown', 'brown', 'brown', 'red', 'purple', 'cyan', 'brown', 'cyan', 'cyan', 'cyan', 'brown', 'brown', 'brown', 'purple', 'brown', 'brown', 'brown', 'brown', 'cyan', 'cyan', 'red', 'green', 'green', 'cyan', 'cyan', 'cyan', 'red', 'cyan', 'cyan', 'brown', 'cyan', 'cyan', 'green', 'cyan', 'orange', 'brown', 'cyan', 'brown', 'brown', 'red', 'black', 'purple', 'brown', 'brown', 'black', 'black', 'brown', 'brown', 'black', 'brown', 'brown', 'brown', 'orange', 'brown', 'brown', 'orange', 'orange', 'green', 'black', 'brown', 'black', 'orange', 'green', 'brown', 'brown', 'brown', 'red', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'green', 'green', 'red', 'cyan', 'brown', 'brown', 'black', 'brown', 'brown', 'black', 'brown', 'purple', 'purple', 'purple', 'purple', 'purple'
)) %>% 
    set("labels_cex", c(0.35)) %>% 
    set("leaves_cex", c(0.45)) %>%
    set("leaves_pch", 19) %>% set("leaves_col", c('red', 'black', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'black', 'brown', 'red', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'black', 'brown', 'brown', 'brown', 'black', 'brown', 'brown', 'red', 'green', 'red', 'brown', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'red', 'red', 'red', 'brown', 'red', 'red', 'brown', 'red', 'red', 'red', 'red', 'green', 'brown', 'brown', 'red', 'red', 'red', 'red', 'red', 'green', 'red', 'brown', 'red', 'red', 'black', 'red', 'black', 'red', 'black', 'red', 'brown', 'brown', 'brown', 'brown', 'red', 'red', 'brown', 'red', 'black', 'red', 'brown', 'brown', 'brown', 'cyan', 'red', 'brown', 'red', 'red', 'red', 'brown', 'red', 'red', 'purple', 'brown', 'red', 'red', 'red', 'red', 'brown', 'brown', 'brown', 'red', 'orange', 'orange', 'orange', 'orange', 'green', 'red', 'brown', 'red', 'green', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'brown', 'brown', 'brown', 'red', 'red', 'red', 'black', 'brown', 'orange', 'brown', 'brown', 'brown', 'brown', 'red', 'red', 'black', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'orange', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'brown', 'red', 'red', 'green', 'red', 'red', 'red', 'green', 'green', 'green', 'green', 'green', 'green', 'green', 'green', 'brown', 'orange', 'green', 'green', 'brown', 'black', 'brown', 'brown', 'green', 'green', 'green', 'green', 'green', 'green', 'green', 'brown', 'green', 'green', 'brown', 'green', 'green', 'green', 'green', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'red', 'orange', 'red', 'brown', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'black', 'black', 'brown', 'brown', 'brown', 'brown', 'brown', 'black', 'brown', 'brown', 'purple', 'black', 'red', 'black', 'green', 'brown', 'brown', 'brown', 'brown', 'cyan', 'brown', 'black', 'cyan', 'orange', 'brown', 'brown', 'green', 'green', 'red', 'red', 'black', 'brown', 'brown', 'red', 'brown', 'green', 'green', 'green', 'brown', 'brown', 'red', 'purple', 'red', 'brown', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'orange', 'purple', 'brown', 'brown', 'purple', 'orange', 'brown', 'cyan', 'brown', 'green', 'brown', 'brown', 'purple', 'purple', 'brown', 'purple', 'brown', 'purple', 'purple', 'purple', 'brown', 'black', 'brown', 'brown', 'brown', 'brown', 'red', 'purple', 'cyan', 'brown', 'cyan', 'cyan', 'cyan', 'brown', 'brown', 'brown', 'purple', 'brown', 'brown', 'brown', 'brown', 'cyan', 'cyan', 'red', 'green', 'green', 'cyan', 'cyan', 'cyan', 'red', 'cyan', 'cyan', 'brown', 'cyan', 'cyan', 'green', 'cyan', 'orange', 'brown', 'cyan', 'brown', 'brown', 'red', 'black', 'purple', 'brown', 'brown', 'black', 'black', 'brown', 'brown', 'black', 'brown', 'brown', 'brown', 'orange', 'brown', 'brown', 'orange', 'orange', 'green', 'black', 'brown', 'black', 'orange', 'green', 'brown', 'brown', 'brown', 'red', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'green', 'green', 'red', 'cyan', 'brown', 'brown', 'black', 'brown', 'brown', 'black', 'brown', 'purple', 'purple', 'purple', 'purple', 'purple'
))


par(mar = rep(0, 4))
circlize_dendrogram(dend, dend_track_height = 0.8)

################################## 
library(StAMPP)	
potato.freq <- stamppConvert("GS_QC_AABB_110g_6575m.csv")
 
# calculate pair-wise Nei’s genetic distances between individuals (FALSE)
potato.ind.d <- stamppNeisD(potato.freq, pop=FALSE)
 
library("factoextra")
obj <- as.dist(potato.ind.d)
 
library(ape)
library(phangorn)
library(seqinr)
library(dendextend)
library(factoextra)
library(circlize)


dend <- obj %>% 
    hclust(method = "ward.D2") %>% as.dendrogram %>%
    set("branches_k_color", value = c("cyan", "blue", "orange"), k=3) %>% set("branches_lwd", 2.2) %>%
    set("labels_cex", c(.9,0.75)) %>% 
    set("labels_colors", c('brown', 'brown', 'green', 'red', 'brown', 'black', 'black', 'red', 'red', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'black', 'red', 'brown', 'black', 'black', 'brown', 'red', 'red', 'cyan', 'red', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'red', 'red', 'brown', 'brown', 'green', 'brown', 'red', 'red', 'brown', 'brown', 'red', 'red', 'red', 'red', 'brown', 'brown', 'black', 'brown', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'red', 'brown', 'brown', 'brown', 'black', 'brown', 'red', 'black', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'purple', 'purple', 'purple', 'purple', 'purple', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'black', 'brown', 'red', 'black', 'purple', 'purple', 'purple', 'black', 'brown', 'black', 'black', 'brown', 'black', 'green', 'green', 'black', 'brown', 'black', 'green', 'purple', 'brown', 'brown', 'yellow', 'brown')) %>% 
    set("labels_cex", c(0.75)) %>% 
    set("leaves_cex", c(1.25)) %>%
    set("leaves_pch", 19) %>% set("leaves_col", c('brown', 'brown', 'green', 'red', 'brown', 'black', 'black', 'red', 'red', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'black', 'red', 'brown', 'black', 'black', 'brown', 'red', 'red', 'cyan', 'red', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'red', 'red', 'brown', 'brown', 'green', 'brown', 'red', 'red', 'brown', 'brown', 'red', 'red', 'red', 'red', 'brown', 'brown', 'black', 'brown', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'red', 'brown', 'brown', 'brown', 'black', 'brown', 'red', 'black', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'purple', 'purple', 'purple', 'purple', 'purple', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'black', 'brown', 'red', 'black', 'purple', 'purple', 'purple', 'black', 'brown', 'black', 'black', 'brown', 'black', 'green', 'green', 'black', 'brown', 'black', 'green', 'purple', 'brown', 'brown', 'yellow', 'brown'))

par(mar = rep(0, 4))
circlize_dendrogram(dend, dend_track_height = 0.8)

########################################### 
import pandas as pd
# convert excel file to .csv file
read_file = pd.read_excel("HC_order_labels3.xlsx",
                          "Blad1")
read_file.to_csv("colors_393m.csv", index=False)

# CSV converted to list with string elements
string1 = "brown,brown,green,red,brown,brown,red,cyan,red,red,red,brown,red,red,red,black,green,red,brown,red,red,red,brown,brown,red,brown,brown,red,red,red,brown,black,orange,brown,brown,red,red,red,red,red,red,red,red,red,brown,brown,red,green,red,red,brown,brown,red,brown,orange,red,red,red,red,red,red,red,red,red,red,red,red,red,red,red,red,red,brown,red,red,green,red,red,red,brown,brown,black,brown,brown,orange,brown,green,red,orange,orange,orange,brown,red,brown,black,brown,brown,black,purple,red,brown,purple,purple,green,green,brown,orange,brown,brown,red,red,red,red,red,brown,red,red,red,green,green,red,brown,brown,brown,red,brown,red,brown,red,red,red,orange,orange,brown,red,brown,brown,brown,red,green,brown,brown,black,red,red,orange,black,red,brown,red,brown,red,red,red,black,red,brown,brown,red,red,red,red,red,red,red,red,brown,orange,brown,brown,brown,brown,orange,brown,orange,purple,orange,purple,purple,purple,orange,brown,brown,brown,brown,green,green,green,green,green,green,green,brown,green,green,green,green,brown,green,green,green,red,brown,brown,brown,red,brown,red,black,black,brown,brown,brown,brown,brown,brown,brown,brown,brown,brown,brown,red,brown,brown,brown,black,brown,red,brown,brown,brown,brown,brown,brown,red,orange,red,brown,brown,brown,brown,brown,orange,brown,brown,brown,purple,purple,purple,brown,brown,green,red,green,purple,brown,cyan,brown,brown,brown,black,black,green,orange,brown,purple,brown,orange,orange,purple,brown,orange,brown,brown,brown,green,brown,black,purple,brown,brown,black,brown,cyan,cyan,brown,cyan,red,green,cyan,cyan,cyan,cyan,cyan,brown,cyan,cyan,green,cyan,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,purple,cyan,brown,brown,brown,brown,brown,brown,brown,brown,brown,red,purple,cyan,brown,cyan,cyan,cyan,brown,brown,brown,brown,brown,brown,brown,brown,red,black,brown,green,brown,green,brown,green,green,brown,green,brown,cyan,brown,orange,black,brown,brown,brown,brown,brown,brown,orange,brown,brown,green,green,green,green,green,brown,brown,green,brown,brown,purple,purple,purple,purple,purple,purple,purple"
print(string1.split(','))	

#output
'brown', 'brown', 'green', 'red', 'brown', 'brown', 'red', 'cyan', 'red', 'red', 'red', 'brown', 'red', 'red', 'red', 'black', 'green', 'red', 'brown', 'red', 'red', 'red', 'brown', 'brown', 'red', 'brown', 'brown', 'red', 'red', 'red', 'brown', 'black', 'orange', 'brown', 'brown', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'brown', 'brown', 'red', 'green', 'red', 'red', 'brown', 'brown', 'red', 'brown', 'orange', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'brown', 'red', 'red', 'green', 'red', 'red', 'red', 'brown', 'brown', 'black', 'brown', 'brown', 'orange', 'brown', 'green', 'red', 'orange', 'orange', 'orange', 'brown', 'red', 'brown', 'black', 'brown', 'brown', 'black', 'purple', 'red', 'brown', 'purple', 'purple', 'green', 'green', 'brown', 'orange', 'brown', 'brown', 'red', 'red', 'red', 'red', 'red', 'brown', 'red', 'red', 'red', 'green', 'green', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'red', 'brown', 'red', 'red', 'red', 'orange', 'orange', 'brown', 'red', 'brown', 'brown', 'brown', 'red', 'green', 'brown', 'brown', 'black', 'red', 'red', 'orange', 'black', 'red', 'brown', 'red', 'brown', 'red', 'red', 'red', 'black', 'red', 'brown', 'brown', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'red', 'brown', 'orange', 'brown', 'brown', 'brown', 'brown', 'orange', 'brown', 'orange', 'purple', 'orange', 'purple', 'purple', 'purple', 'orange', 'brown', 'brown', 'brown', 'brown', 'green', 'green', 'green', 'green', 'green', 'green', 'green', 'brown', 'green', 'green', 'green', 'green', 'brown', 'green', 'green', 'green', 'red', 'brown', 'brown', 'brown', 'red', 'brown', 'red', 'black', 'black', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'brown', 'brown', 'brown', 'black', 'brown', 'red', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'orange', 'red', 'brown', 'brown', 'brown', 'brown', 'brown', 'orange', 'brown', 'brown', 'brown', 'purple', 'purple', 'purple', 'brown', 'brown', 'green', 'red', 'green', 'purple', 'brown', 'cyan', 'brown', 'brown', 'brown', 'black', 'black', 'green', 'orange', 'brown', 'purple', 'brown', 'orange', 'orange', 'purple', 'brown', 'orange', 'brown', 'brown', 'brown', 'green', 'brown', 'black', 'purple', 'brown', 'brown', 'black', 'brown', 'cyan', 'cyan', 'brown', 'cyan', 'red', 'green', 'cyan', 'cyan', 'cyan', 'cyan', 'cyan', 'brown', 'cyan', 'cyan', 'green', 'cyan', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'cyan', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'purple', 'cyan', 'brown', 'cyan', 'cyan', 'cyan', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'red', 'black', 'brown', 'green', 'brown', 'green', 'brown', 'green', 'green', 'brown', 'green', 'brown', 'cyan', 'brown', 'orange', 'black', 'brown', 'brown', 'brown', 'brown', 'brown', 'brown', 'orange', 'brown', 'brown', 'green', 'green', 'green', 'green', 'green', 'brown', 'brown', 'green', 'brown', 'brown', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple', 'purple'

########################################## 



















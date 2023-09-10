###--------------------###
bash Anaconda3-2021.05-Linux-x86_64.sh

conda install -c dranew shapeit

1. download shapeit.v2.r904.glibcv2.12.linux.tar.gz
2. open file(shapeit.v2.904.3.10.0-693.11.6.el7.x86_64)
3. move four files (gwas.bed, gwas.bim, gwas.fam and genetic map.txt) in the example folder to the anaconda3/bin
4. move a file shapeit to the anaconda3/bin
5. open terminal in the anaconda3/bin

shapeit --input-gen chip_chr11.gen chip_potato_chr11.sample \
        -M genetic_map_chr11.txt \
        -O chip_chr11.phased \
        --force \
        --window 0.5 \
        --burn 80 \
        --prune 60 \
        --main 100
###---------------------###


library(rehh)

# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 1, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr01.map", 
                   chr.name = "chr01",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_1 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 2, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr02.map", 
                   chr.name = "chr02",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_2 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 3, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr03.map", 
                   chr.name = "chr03",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_3 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 4, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr04.map", 
                   chr.name = "chr04",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_4 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 5, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr05.map", 
                   chr.name = "chr05",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_5 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 6, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr06.map", 
                   chr.name = "chr06",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_6 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 7, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr07.map", 
                   chr.name = "chr07",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_7 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 8, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr08.map", 
                   chr.name = "chr08",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_8 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 9, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr09.map", 
                   chr.name = "chr09",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_9 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 10, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr10.map", 
                   chr.name = "chr10",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_10 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 11, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr11.map", 
                   chr.name = "chr11",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_11 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/SK/hap_chr_", 12, ".sk", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/SK/chr12.map", 
                   chr.name = "chr12",
                   allele_coding = "01")

# perform scan on a single chromosome (calculate iHH values)
scan_12 <- scan_hh(hh)


wgscan.sk <- rbind(scan_1, scan_2, scan_3, scan_4, scan_5, scan_6, scan_7, scan_8, scan_9, scan_10, scan_11, scan_12)
# calculate genome-wide iHS values
wgscan.ihs.sk <- ihh2ihs(wgscan.sk)	

head(wgscan.ihs.sk$ihs)
head(wgscan.ihs.sk$frequency.class)
library(rehh)


# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 1, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr01.map", 
                   chr.name = "chr01",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_1 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 2, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr02.map", 
                   chr.name = "chr02",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_2 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 3, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr03.map", 
                   chr.name = "chr03",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_3 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 4, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr04.map", 
                   chr.name = "chr04",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_4 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 5, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr05.map", 
                   chr.name = "chr05",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_5 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 6, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr06.map", 
                   chr.name = "chr06",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_6 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 7, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr07.map", 
                   chr.name = "chr07",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_7 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 8, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr08.map", 
                   chr.name = "chr08",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_8 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 9, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr09.map", 
                   chr.name = "chr09",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_9 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 10, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr10.map", 
                   chr.name = "chr10",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_10 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 11, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr11.map", 
                   chr.name = "chr11",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_11 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/JP/hap_chr_", 12, ".jp", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/JP/chr12.map", 
                   chr.name = "chr12",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_12 <- scan_hh(hh)


wgscan.jp <- rbind(scan_1, scan_2, scan_3, scan_4, scan_5, scan_6, scan_7, scan_8, scan_9, scan_10, scan_11, scan_12)
# calculate genome-wide iHS values
wgscan.ihs.jp <- ihh2ihs(wgscan.jp)	

head(wgscan.ihs.jp$ihs)
head(wgscan.ihs.jp$frequency.class)
rsb.sk_jp <- ines2rsb(scan_pop1 = wgscan.sk,
                        scan_pop2 = wgscan.jp,
                        popname1 = "SK",
                        popname2 = "JP")
Scan of pop1 contains 4004 markers.
Scan of pop2 contains 4004 markers.
Merged data contains 4004 markers.
head(rsb.sk_jp)
write.csv(rsb.sk_jp, file = "D:/WORK/NEW Project/25-SHAPEIT2/rsb.sk_jp.csv")
# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 1, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr01.map", 
                   chr.name = "chr01",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_1 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 2, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr02.map", 
                   chr.name = "chr02",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_2 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 3, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr03.map", 
                   chr.name = "chr03",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_3 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 4, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr04.map", 
                   chr.name = "chr04",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_4 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 5, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr05.map", 
                   chr.name = "chr05",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_5 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 6, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr06.map", 
                   chr.name = "chr06",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_6 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 7, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr07.map", 
                   chr.name = "chr07",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_7 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 8, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr08.map", 
                   chr.name = "chr08",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_8 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 9, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr09.map", 
                   chr.name = "chr09",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_9 <- scan_hh(hh)


# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 10, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr10.map", 
                   chr.name = "chr10",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_10 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 11, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr11.map", 
                   chr.name = "chr11",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_11 <- scan_hh(hh)

# haplotype file name for each chromosome
hap_file = paste("D:/US/hap_chr_", 12, ".us", sep = "")
# create internal representation
hh <- data2haplohh(hap_file = hap_file,
                   map_file = "D:/US/chr12.map", 
                   chr.name = "chr12",
                   allele_coding = "01")
				   
# perform scan on a single chromosome (calculate iHH values)
scan_12 <- scan_hh(hh)


wgscan.us <- rbind(scan_1, scan_2, scan_3, scan_4, scan_5, scan_6, scan_7, scan_8, scan_9, scan_10, scan_11, scan_12)
# calculate genome-wide iHS values
wgscan.ihs.us <- ihh2ihs(wgscan.us)	

head(wgscan.ihs.us$ihs)
head(wgscan.ihs.us$frequency.class)
rsb.sk_us <- ines2rsb(scan_pop1 = wgscan.sk,
                        scan_pop2 = wgscan.us,
                        popname1 = "SK",
                        popname2 = "US")
head(rsb.sk_us)
write.csv(rsb.sk_us, file = "D:/WORK/NEW Project/25-SHAPEIT2/rsb.sk_us.csv")
xpehh.jp_us <- ies2xpehh(scan_pop1 =  wgscan.jp,
                         scan_pop2 =  wgscan.us,
						 popname1 = "JP",
						 popname2 = "US")

head(xpehh.jp_us)
write.csv(xpehh.jp_us, file = "D:/WORK/NEW Project/25-SHAPEIT2/xpehh.jp_us.csv")
# Rsb vs. XP-EHH comparison
plot(rsb.sk_us[, "RSB_SK_US"],
     xpehh.sk_us[, "XPEHH_SK_US"],
     xlab = "Rsb",
     ylab = "XP-EHH",
     pch = ".",
     xlim = c(-7.5, 7.5),
     ylim = c(-7.5, 7.5))
# add red circle for marker with name "solcap_snp_c2_26799"
points(rsb.sk_us["solcap_snp_c2_26799", "RSB_SK_US"],
       xpehh.sk_us["solcap_snp_c2_26799", "XPEHH_SK_US"],
       col = "red")
# add dashed diagonal
abline(a = 0, b = 1, lty = 2)
 
distribplot(wgscan.ihs.sk$ihs$IHS, xlab = "iHS")
    
distribplot(wgscan.ihs.sk$ihs$IHS, 
            xlab = "iHS", 
            qqplot = TRUE)
±
# Rsb vs. XP-EHH comparison
plot(rsb.jp_us[, "RSB_JP_US"],
     xpehh.jp_us[, "XPEHH_JP_US"],
     xlab = "Rsb",
     ylab = "XP-EHH",
     pch = ".",
     xlim = c(-7.5, 7.5),
     ylim = c(-7.5, 7.5))
# add red circle for marker with name "solcap_snp_c2_26799"
points(rsb.jp_us["solcap_snp_c2_26799", "RSB_JP_US"],
       xpehh.jp_us["solcap_snp_c2_26799", "XPEHH_JP_US"],
       col = "red")
# add dashed diagonal
abline(a = 0, b = 1, lty = 2)
 

manhattanplot(xpehh.jp_us,
              threshold = 3,
			  #threshold = -3,
              main = "XPEHH Japanese potatoes_US potatoes")

 
manhattanplot(xpehh.jp_us,
              threshold = c(3,-3),
              main = "XPEHH Japanese potatoes_US potatoes")
 
manhattanplot(xpehh.jp_us,
              pval = TRUE,
              threshold = 3,
              main = "p-value of XPEHH (Japanese potatoes_US potatoes)")

 



cr.cgu <- calc_candidate_regions(xpehh.jp_us,
                                 threshold = 3,
                                 pval = TRUE,
                                 window_size = 1E6,
                                 overlap = 1E5,
                                 min_n_extr_mrk = 2)

# re-define colors
palette(c("red", "green"))
manhattanplot(xpehh.jp_us, 
              pval = TRUE,
              threshold = 3, 
              chr.name = c("chr02", "chr04", "chr07", "chr08"), 
              main = "XPEHH (Japanese potatoes_US potatoes)", 
              cr = cr.cgu,
              mrk = "loc_3745",
              inset = 1E+7,
              resolution = c(200000, 0.05))
# Rasterization reduced 1235 data points to 998 .
 
# set back to default colors
palette("default")
 


xpehh.sk_jp <- ies2xpehh(scan_pop1 =  wgscan.sk,
                         scan_pop2 =  wgscan.jp,
                         popname1 = "SK",
                         popname2 = "JP")

manhattanplot(xpehh.sk_jp,
              threshold = c(3,-3),
              main = "XPEHH Korean potatoes_Japanese potatoes")
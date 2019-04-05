#Example pipeline
rm(list=ls())
gc()



#Set wd to folder containing all data
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/GISTIC/gdac.broadinstitute.org_ACC-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0")

focal_input <- read.table("focal_input.seg.txt",header=TRUE)
cutoff <- read.table("sample_cutoffs.txt",header=TRUE)
gene_level1 <- read.table("all_data_by_genes.txt",fill=TRUE,header=TRUE)

gene_tresholded <- read.table("all_thresholded.by_genes.txt",fill=TRUE,header=TRUE)



gene_level <- gene_level1[, 6:ncol(gene_level1)]


gene_level[gene_level < -1.054] <- -2
gene_level[gene_level > 0.336] <- 2
gene_level[gene_level > 0.1 & gene_level <= 0.336 ] <- 1
gene_level[gene_level < -0.1 & gene_level >= -1.054 ] <- -1
gene_level[gene_level >= -0.1 & gene_level <= 0.1 ] <- 0

gene_level1[, 6:ncol(gene_level1)] <- gene_level


test <- head(gene_tresholded[, 6:ncol(gene_tresholded)]) - head(gene_level)



focal_input <- focal_input %>% filter(Sample == "TCGA-OR-A5J1-01A-11D-A29H-01")

CN_vector <- focal_input$Seg.CN


CN_vector[CN_vector <= -1.054] <- -2
CN_vector[CN_vector >= 0.336] <- 2
CN_vector[CN_vector >= 0.1 & CN_vector < 0.336 ] <- 1
CN_vector[CN_vector <= -0.1 & CN_vector > -1.054 ] <- -1
CN_vector[CN_vector > -0.1 & CN_vector < 0.1 ] <- 0

focal_input$tresholded <- CN_vector

#what are the genes in these 

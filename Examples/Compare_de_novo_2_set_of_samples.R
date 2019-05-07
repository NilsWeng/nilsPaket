#Compare two set of samples , extract de novo-signatures for both of them

rm(list=ls())
gc()
#Library
library(nilsPaket)
library(dplyr)
library(MutationalPatterns)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(gridExtra)
library(grid)


#Set wd to folder containing all data
main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))
#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"


#Vector with TCGA-barcodes
#set1 <- "BLA"
#set2 <- "BLI"





# Set1 --------------------------------------------------------------------

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output")

file_name <- "SD_random3_mut_mat.rda"
load(file_name)
random_mut_mat <- mutational_matrix
colnames(random_mut_mat) <- c(1:ncol(random_mut_mat))
#colnames(random_mut_mat) <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(random_mut_mat))


#Find rank optimisation
random_mut_mat <- random_mut_mat + 0.0001
estimate <- nmf(random_mut_mat, rank=2:15, method="brunet", nrun=10, seed=123456)
plot(estimate)


#Extract de novo signatures

number_of_signatures <- 6
de_novo_random <- extract_signatures(random_mut_mat, rank = number_of_signatures, nrun = 10)

colnames(de_novo_random$signatures) <- c(1:number_of_signatures)
rownames(de_novo_random$contribution) <- c(1:number_of_signatures)




plot_96_profile(de_novo_random$signatures, condensed = TRUE)


colnames(de_novo_random$signatures) <- c(paste(rep("De-novo",number_of_signatures),c(1:number_of_signatures),sep=" "))
rownames(de_novo_random$contribution) <- c(paste(rep("De-novo",number_of_signatures),c(1:number_of_signatures),sep=" "))



#Plot
pc1 <- plot_contribution(de_novo_random$contribution, de_novo_random$signature,mode = "relative")
#Visualize the contribution of the signatures in absolute number of mutations:
pc2 <- plot_contribution(de_novo_random$contribution, de_novo_random$signature,mode = "absolute")
#Combine the two plots:
grid.arrange(pc1, pc2)


#Cosine similiarity between de-novo and extracted
cos_similarity_signatures <- cos_sim_matrix(de_novo_random$signatures, cosmic_signatures)

hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]


plot_cosine_heatmap(cos_similarity_signatures,col_order = cosmic_order,cluster_rows = FALSE)





# Set2 --------------------------------------------------------------------


setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output")

file_name <- "SD_3_mut_mat.rda"
load(file_name)
hyper_mutational_matrix <- mutational_matrix

#colnames(hyper_mutational_matrix) <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(hyper_mutational_matrix))
colnames(hyper_mutational_matrix) <- c(1:ncol(hyper_mutational_matrix))

#Find rank optimisation
hyper_mutational_matrix <- hyper_mutational_matrix + 0.0001
estimate <- nmf(hyper_mutational_matrix, rank=2:25, method="brunet", nrun=10, seed=123456)
plot(estimate)


#Extract de novo signatures

number_of_signatures <- 6
de_novo_hyper <- extract_signatures(hyper_mutational_matrix, rank = number_of_signatures, nrun = 10)

colnames(de_novo_hyper$signatures) <- c(1:number_of_signatures)
rownames(de_novo_hyper$contribution) <- c(1:number_of_signatures)




plot_96_profile(de_novo_hyper$signatures, condensed = TRUE)


colnames(de_novo_hyper$signatures) <- c(paste(rep("De-novo",number_of_signatures),c(1:number_of_signatures),sep=" "))
rownames(de_novo_hyper$contribution) <- c(paste(rep("De-novo",number_of_signatures),c(1:number_of_signatures),sep=" "))



#Plot
pc1 <- plot_contribution(de_novo_hyper$contribution, de_novo_hyper$signature,mode = "relative")
#Visualize the contribution of the signatures in absolute number of mutations:
pc2 <- plot_contribution(de_novo_hyper$contribution, de_novo_hyper$signature,mode = "absolute")
#Combine the two plots:
grid.arrange(pc1, pc2)


#Cosine similiarity between de-novo and extracted
cos_similarity_signatures <- cos_sim_matrix(de_novo_hyper$signatures, cosmic_signatures)

hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]


plot_cosine_heatmap(cos_similarity_signatures,col_order = cosmic_order,cluster_rows = FALSE)





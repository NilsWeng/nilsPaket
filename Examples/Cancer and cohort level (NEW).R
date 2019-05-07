#Cancer and cohort level 
rm(list=ls())
gc()
#Library
library(nilsPaket)
library(dplyr)
library(MutationalPatterns)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)



#Set wd to folder containing all data
main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)


#Read cosmic signatures
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))
#cosmic_signatures <- as.matrix(read.table("cosmic_signatures.txt",header=TRUE))


#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output")
#mutational_matrix <- read.table("COMPLETE_MUT_MAT.txt",header=TRUE)


setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output/NEW_COHORT")
mut_mat_files <- list.files()
mut_mat_files <- mut_mat_files[grep("_mut_matrix.txt",mut_mat_files)]




#mut_mat_files <- mut_mat_files[1]


contribution_DF <- data.frame()

for (cancer_file in mut_mat_files){
  
  cancer <- gsub("_..*$","",cancer_file)
  
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output/NEW_COHORT")
 
  cancer_mut_mat <- read.table(cancer_file,header=TRUE)
  colnames(cancer_mut_mat) <- gsub("\\.","-",colnames(cancer_mut_mat))
  
  #Normalize each column by its colsums 
  #normalized_mut_mat <- sweep(cancer_mut_mat,2,colSums(cancer_mut_mat),`/`)

  
  #Summarize into one big cohort normalized mutational profile
  #cohort_normalized_mut_mat <- rowSums(normalized_mut_mat)
  
  
  #Alternativley remove high mut samples from samples then rowsums
  
  #test <- colSums(cancer_mut_mat)
  test <- cancer_mut_mat[which(colSums(cancer_mut_mat) >  2000 )]
  #test <- cancer_mut_mat
  test <- rowSums(test)
  test <- as.matrix(test)
  
  fit_to_cosmic <- fit_to_signatures(test,cosmic_signatures)
  
  contribution <- fit_to_cosmic$contribution
  contribution <- contribution[order(contribution,decreasing = TRUE),]
  
  top_x_signatures <- head(contribution,n=10)
  
  append_df <- data.frame("Contribution"=as.vector(top_x_signatures),"Signature"=names(top_x_signatures))
  
  append_df$cancer <- rep(cancer,nrow(append_df))
  
  contribution_DF <- rbind(contribution_DF,append_df)
  
  
  
  

  
}

contribution_DF <- contribution_DF %>% filter(Contribution > 0)

test <- dplyr::count(contribution_DF,Signature)





# Extract de-novo signatures for each cancer  -----------------------------

rm(list=ls())
gc()
#Library
library(nilsPaket)
library(dplyr)
library(MutationalPatterns)
library(BSgenome)
library(gridExtra)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library("NMF")


main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output/NEW_COHORT")
mut_mat_files <- list.files()
mut_mat_files <- mut_mat_files[grep("_mut_matrix.txt",mut_mat_files)]



cohort_mut_mat <- list()
#to_keep <- c(1,2,3,4,5,6,7,8)

for (file_name in mut_mat_files){
  
  
  
  cancer <- gsub("_..*$","",file_name)
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output/NEW_COHORT")
  cancer_mut_mat <- read.table(file_name,header=TRUE)
  colnames(cancer_mut_mat) <- gsub("\\.","-",colnames(cancer_mut_mat))

  
  #Filter away high_mut samples
  if(TRUE){
    
    number_mut_vec <- as.numeric(colSums(cancer_mut_mat))
    m <- median(number_mut_vec)
    s <- sd(number_mut_vec)
    
    to_keep <- which(number_mut_vec < (m+3*s))
    to_keep <- names(cancer_mut_mat)[to_keep]
    
    cancer_mut_mat <- cancer_mut_mat[,(names(cancer_mut_mat) %in% to_keep)]
    
    
    #to_discard <- which(number_mut_vec < (m+3*s))
    #to_discard <- names(cancer_mut_mat)[to_discard]
    
    
    #cancer_mut_mat <- cancer_mut_mat[, to_keep]

    #cancer_mut_mat <- cancer_mut_mat[,!(names(cancer_mut_mat) %in% to_discard)]
    
  }
  
  
  
  
  if (length(to_keep) > 1){
    
    cohort_mut_mat[[cancer]] <- rowSums(cancer_mut_mat)
    
    
  }else {#If only one sample is in to_keep
    
    cohort_mut_mat[[cancer]] <- cancer_mut_mat
    print("HELLO")
  }
    
  
  
  
}


cohort_mut_mat <- do.call(rbind, cohort_mut_mat)
cohort_mut_mat <- t(cohort_mut_mat)




#setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Pictures")
#pdf("De-novo.pdf",width=16, height=10)

#Find rank optimisation
#cohort_mut_mat <- cohort_mut_mat + 0.0001
#estimate <- nmf(cohort_mut_mat, rank=2:15, method="brunet", nrun=10, seed=123456)
#plot(estimate)


#Extract de novo signatures

number_of_signatures <- 5
de_novo <- extract_signatures(cohort_mut_mat, rank = number_of_signatures, nrun = 10)

colnames(de_novo$signatures) <- c(1:number_of_signatures)
rownames(de_novo$contribution) <- c(1:number_of_signatures)




plot_96_profile(de_novo$signatures, condensed = TRUE)


colnames(de_novo$signatures) <- c(paste(rep("De-novo",number_of_signatures),c(1:number_of_signatures),sep=" "))
rownames(de_novo$contribution) <- c(paste(rep("De-novo",number_of_signatures),c(1:number_of_signatures),sep=" "))



#Plot
pc1 <- plot_contribution(de_novo$contribution, de_novo$signature,mode = "relative")
#Visualize the contribution of the signatures in absolute number of mutations:
pc2 <- plot_contribution(de_novo$contribution, de_novo$signature,mode = "absolute")
#Combine the two plots:
grid.arrange(pc1, pc2)


#Cosine similiarity between de-novo and extracted
cos_similarity_signatures <- cos_sim_matrix(de_novo$signatures, cosmic_signatures)

hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]


plot_cosine_heatmap(cos_similarity_signatures,col_order = cosmic_order,cluster_rows = FALSE)


#dev.off()
#de_novo$signatures <-  de_novo$signatures[c(4,5,1,3,2)]
#de_novo$contribution

test1 <- plot_96_profile(de_novo_high$signatures, condensed = TRUE)
test2 <- plot_96_profile(de_novo$signatures, condensed = TRUE)
test3 <- plot_contribution(de_novo_high$contribution, de_novo_high$signature,mode = "absolute")
test4 <- plot_contribution(de_novo$contribution, de_novo$signature,mode = "absolute")
#test5 <- 
#test6 <- 

grid.arrange(test1,test2,test3,test4,ncol=2,nrow=2)


# Extract for all samples -------------------------------------------------

rm(list=ls())
gc()
#Library
library(nilsPaket)
library(dplyr)
library(MutationalPatterns)
library(BSgenome)
library(gridExtra)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library("NMF")


main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/output/NEW_COHORT")
mut_mat_files <- list.files()
mut_mat_files <- mut_mat_files[grep("_mut_matrix.txt",mut_mat_files)]



datalist = lapply(mut_mat_files,function(x){read.table(file=x,header=TRUE)})

total_mut_mat <- do.call(cbind, datalist)
colnames(total_mut_mat) <- gsub("\\.","-",colnames(total_mut_mat))

total_mut_mat <- as.matrix(total_mut_mat)

cohort_mut_mat <- total_mut_mat



number_of_signatures <- 20
de_novo <- extract_signatures(cohort_mut_mat, rank = number_of_signatures, nrun = 10)


setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data")
save(de_novo,file="de_novo_all.rda")



colnames(de_novo$signatures) <- c(1:number_of_signatures)
rownames(de_novo$contribution) <- c(1:number_of_signatures)




plot_96_profile(de_novo$signatures, condensed = TRUE)


colnames(de_novo$signatures) <- c(paste(rep("De-novo",number_of_signatures),c(1:number_of_signatures),sep=" "))
rownames(de_novo$contribution) <- c(paste(rep("De-novo",number_of_signatures),c(1:number_of_signatures),sep=" "))



#Plot
pc1 <- plot_contribution(de_novo$contribution, de_novo$signature,mode = "relative")
#Visualize the contribution of the signatures in absolute number of mutations:
pc2 <- plot_contribution(de_novo$contribution, de_novo$signature,mode = "absolute")
#Combine the two plots:
grid.arrange(pc1, pc2)


#Cosine similiarity between de-novo and extracted
cos_similarity_signatures <- cos_sim_matrix(de_novo$signatures, cosmic_signatures)

hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
# store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]


plot_cosine_heatmap(cos_similarity_signatures,col_order = cosmic_order,cluster_rows = TRUE)

dev.off()




#total_mut_mat <-  t(total_mut_mat)  #### Dimensioner?






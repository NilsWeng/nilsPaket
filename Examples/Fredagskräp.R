#Example pipeline
rm(list=ls())
gc()
#Library
library(nilsPaket)
library(dplyr)
library(MutationalPatterns)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19", character.only = TRUE)
library(ggplot2)
library(reshape2)



#Set wd to folder containing all data
main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)


#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))


#Get all mutational matrix
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/generated_data")
mut_files <- list.files(pattern = "mut_mat",recursive = TRUE)
mut_files <- mut_files[3]

#for (file_name in mut_files){
  file_name <- mut_files
  
  #load mutational matrix
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Data/MC3/generated_data")
  load(file_name)
  cancer <- gsub("/..*","",file_name)
  
 
  
  
  #Saved in wierd format before need to restore TCGA-name
  setwd(main_wd)
  load("MC3.rda")
  sample_names <- unique(MC3$Tumor_Sample_Barcode)
  mat_names <- colnames(mutational_matrix)
  colnames(mutational_matrix) <- grep(paste(mat_names, collapse="|"), sample_names,value = TRUE)
  
  

  
  #Create a cosine similarity matrix
  cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix, cosmic_signatures)
  
  
  #Fit mut_matrix to cosmic
  fit_to_cosmic <- fit_to_signatures(mutational_matrix,cosmic_signatures)
  

  #Try cluster data before plotting
  sample_order <- hclust(dist(t(fit_to_cosmic$contribution)),method="complete")$order
  sample_order <- colnames(fit_to_cosmic$contribution)[sample_order]
  
  #
  fit_to_cosmic$contribution <- fit_to_cosmic$contribution[, match(sample_order,colnames(fit_to_cosmic$contribution))]
  fit_to_cosmic$reconstructed <- fit_to_cosmic$reconstructed[, match(sample_order,colnames(fit_to_cosmic$reconstructed))]
  mutational_matrix <- mutational_matrix[, match(sample_order,colnames(mutational_matrix))]
  
  
  
  
  #Apply some kind of filter to only get main contributors
  select <- which(rowSums(fit_to_cosmic$contribution)/sum(rowSums(fit_to_cosmic$contribution)) > 0.02)
  #Subselection of only high contributors , cluster on this.
  #test <- fit_to_cosmic$contribution
  #test <- test[rownames(test) %in% sample_order ,]
  
  #Plot contribution
  plot_contribution(fit_to_cosmic$contribution[select ,],cosmic_signatures[,select],coord_flip= FALSE,mode="absolute")
  

 


  contribution <- fit_to_cosmic$contribution[select ,]
  Signature <- cosmic_signatures[,select]
  
  m_contribution = melt(contribution)
  colnames(m_contribution) = c("Signature", "Sample", "Contribution")
  
  plot = ggplot(m_contribution,
                aes(x = factor(Sample),
                    y = Contribution,
                    fill = factor(Signature),
                    order = Sample)) +
    geom_bar(position = "fill", stat="identity", colour="black")  +
    # ylabel
    labs(x = "", y = "Relative contribution") +
    # white background
    theme_bw() +
    # no gridlines
    theme(panel.grid.minor.x=element_blank(),
          panel.grid.major.x=element_blank()) +
    theme(panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          axis.text.y=element_text(vjust = 0.5,size=5))
  

plot = plot + scale_fill_discrete(name="Signature")
plot = plot + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample))))
plot  

  
  

  
#How well is the original signatures reconstructed
cos_sim_reco <- cos_sim_matrix(mutational_matrix, fit_to_cosmic$reconstructed)
cos_sim_reco <- as.data.frame(diag(cos_sim_reco))
colnames(cos_sim_reco) = "cos_sim"
cos_sim_reco$sample = row.names(cos_sim_reco)

ggplot(cos_sim_reco, aes(y=cos_sim, x=sample)) +
   geom_bar(stat="identity", fill = "skyblue4") +
   coord_cartesian(ylim=c(0.6, 1)) +
   coord_flip(ylim=c(0.6,1)) +
   ylab("Cosine similarity\n original VS reconstructed") +
   xlab("") +
   # Reverse order of the samples such that first is up
   # xlim(rev(levels(factor(cos_sim_ori_rec$sample)))) +
   theme_bw() +
   theme(panel.grid.minor.y=element_blank(), panel.grid.major.y=element_blank(),
         axis.text.y=element_text(vjust = 0.5,size=5))+ 
   # Add cut.off line
   geom_hline(aes(yintercept=.95))

  
  
  
  
  
  
  
  
  
  












  
  
  
  #Cluster cosmic signatures
  hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
  #store signatures in new order
  cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]
  #plot cosine_heatmap
  print(plot_cosine_heatmap(cos_sim_samples_cosmic, col_order = cosmic_order, cluster_rows = TRUE))
  
  
  
  #Cluster samples based on cosine similarity
  
  
  
  # N <- 5
  
  sample_cluster <- hclust(dist(cos_sim_samples_cosmic,method="euclidean"),method="complete")
  
  N <- length(unique(cutreeDynamicTree(sample_cluster,maxTreeHeight = 16,deepSplit = FALSE)))
  
  if (N < 3 ){
    
    N <- 3
  }
  
  print(plot_cluster(sample_cluster,N,2.2))
  
  
  
  
  
  #plot cluster in heatmap
  setwd(main_wd)
  ClusterDF <- plot_cluster_in_cosine(sample_cluster,cos_sim_samples_cosmic,N)
  ClusterDF$sample <- as.character(ClusterDF$sample)
  sampleDF$sample <- as.character(sampleDF$sample)
  #ClusterDF$sample <- as.numeric(ClusterDF$sample) 
  ClusterDF<-left_join(ClusterDF,sampleDF,by="sample")
  ClusterDF$TSS.Code <- gsub("TCGA-","",gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",ClusterDF$TCGA_code))
  TSS2Study <- read.table("TSS2Studyabb.txt",header=TRUE)
  ClusterDF <- left_join(ClusterDF,TSS2Study,by="TSS.Code")
  ClusterDF <- ClusterDF %>% select("sample","cluster","TCGA_code","Study.Abbreviation")
  
  
  load("MC3.rda")
  mutDF <- count(MC3,Tumor_Sample_Barcode)
  colnames(mutDF) <- c("TCGA_code","mutations")
  ClusterDF <- left_join(ClusterDF,mutDF,by="TCGA_code")
  
  
  
  rm(TSS2Study,sampleDF,MC3)
  
  
  #Find signatures in cluster
  #Make sure that rownames of cos_sim matches ClusterDF$samples
  #get_signature(ClusterDF,cos_sim_samples_cosmic,0.6)
  
  
  # CNV_SNV plot ------------------------------------------------------------
  
  
  setwd(main_wd)
  
  samples <- ClusterDF %>% select(TCGA_code)
  samples <- as.character(samples$TCGA_code)
  samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples)
  
  
  genes <- read.table("gen_lista.csv",header=TRUE,sep=";")
  genes <- as.character(genes$Gen)
  
  
  print(plot_CNV_SNV(samples,genes,"CNV/SNV plot",plot_id=TRUE,cluster_names = TRUE))
  
  
  
  # mRNA Exp ----------------------------------------------------------------
  
  mRNA_DF <- mRNA_exp(genes)
  
  exp <- get_exp_as_matrix(samples,genes,mRNA_DF)
  exp <- t(exp)
  
  
  
  #For getting cluster lines
  test <- ClusterDF 
  test$TCGA_code <-  gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",ClusterDF$TCGA_code)
  test <- test[which(test$TCGA_code %in% colnames(exp)) ,]
  hline_pos <- c()
  for (cluster1 in unique(test$cluster)){
    
    
    pos <- tail(which(test$cluster == cluster1),n=1)
    pos <- pos + 0.5
    
    hline_pos <- c(hline_pos,pos)
    
    
  }
  hline_pos <-head(hline_pos,-1)
  
  
  
  #CORRECT????? !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  colnames(exp) <- ClusterDF[which(gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",ClusterDF$TCGA_code) %in% colnames(exp)) , 1]
  colnames(exp) <- paste(colnames(exp),"-",sep="")
  
  
  
  
  
  #In order to not get grey values
  exp[exp >2] <- 2
  
  
  exp_melt <- melt(exp)
  
  p <- ggplot(data = exp_melt, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + labs(x="",y="",title="mRNA expression(exp/avg)")+
    #scale_fill_gradient(low="green", mid="lightblue", high="red",limits
    scale_fill_gradientn(colours=c("red","lightblue","green"), limits=c(0,2))+
    geom_hline(yintercept=hline_pos)+
    
    theme(
      axis.text.x=element_text(colour="black",angle=90, hjust=1,vjust = 0.5,size=10),
      axis.text.y=element_text(vjust = 0.5,size=6)
    )
  
  print(p)
  dev.off()
  print(paste("done with ",cancer,sep=""))
  
}


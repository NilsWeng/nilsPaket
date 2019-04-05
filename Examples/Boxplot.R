rm(list=ls())
gc()

library(data.table)
library(dplyr)
library(nilsPaket)
library(ggplot2)
library(MutationalPatterns)        
  
main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)


MC3_DF <- fread('mc3.v0.2.8.PUBLIC.maf.gz')
#save(MC3_DF,file="MC3_DF_tot.rda")
#load("MC3_DF_tot.rda")

Sample_count <- MC3_DF %>% select(Tumor_Sample_Barcode) %>% count(Tumor_Sample_Barcode)
Sample_count$cancer <- get_cancer_type(gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",Sample_count$Tumor_Sample_Barcode))
Cancer_count  <-  Sample_count %>% group_by(cancer) %>% summarise(sum(n))



#Boxplot 
box_p <- ggplot(Sample_count, aes(x=cancer, y=n)) +
  geom_boxplot() +
  stat_boxplot(geom ='errorbar')+
  #theme_bw()+
  labs(title="",x="Cancer", y = "Number of mutations/sample")+
  #geom_boxplot(coef = 6) +
  geom_boxplot(outlier.colour=NA,outlier.size = 0,outlier.alpha = 0) +
  coord_cartesian(ylim = c(0, 3500))+
  theme_minimal()+
  theme(axis.text.x=element_text(vjust = 0.5,size=8,angle = 90))

  
box_p

#In order to remove outliers, since outlier.colour=NA doesnt work
library(plotly)
box_p <- plotly_build(box_p)

for(i in 1:length(box_p$x$data)) {
  box_p$x$data[[i]]$marker$opacity = 0

  }

box_p







# Cosmic contribution -----------------------------------------------------
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/mut_matrix")
load("mut_matrix_cohort.rda")
main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)
load("Sample_count.rda")



#normalize mutational matrix after number of samples in each cohort
samples_in_cancer <- Sample_count %>% select("cancer") %>% dplyr::count(cancer)

for (cancer_type in samples_in_cancer$cancer){
  
  column <- which(colnames(mutational_matrix) %in% cancer_type)
  number_of_samples <- as.integer(samples_in_cancer %>% filter(cancer == cancer_type) %>% select(n))
  mutational_matrix[, column] <- mutational_matrix[, column] / number_of_samples
}
  



#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))

#Create a cosine similarity matrix
cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix, cosmic_signatures)

#Fit mut_matrix to cosmic
fit_to_cosmic <- fit_to_signatures(mutational_matrix,cosmic_signatures)




#Filter signatures to plot
select <- which(rowSums(fit_to_cosmic$contribution)/sum(rowSums(fit_to_cosmic$contribution)) > 0.01)

#Plot
plot_contribution(fit_to_cosmic$contribution[select ,],cosmic_signatures[,select],coord_flip= TRUE,mode="absolute")
plot_cosine_heatmap(cos_sim_samples_cosmic,cluster_rows = TRUE)




#Normalize contribution and plot
test <- fit_to_cosmic$contribution
test <- sweep(test,2,colSums(test),"/")

#Cluster based on contribution similiarity
col_order <- hclust(dist(t(test)),method="complete")$order
col_order <- colnames(test)[col_order]

test <- test[, match(col_order,colnames(test))]
select <- which(rowSums(test) > 1)
plot_contribution(test[select ,],cosmic_signatures[, select],coord_flip= TRUE,mode="absolute")





#Find top X signatures in each cancertype 

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/mut_matrix")
load("mut_matrix_cohort.rda")
setwd(main_wd)

samples_in_cancer <- Sample_count %>% select("cancer") %>% dplyr::count(cancer)

for (cancer_type in samples_in_cancer$cancer){
  
  column <- which(colnames(mutational_matrix) %in% cancer_type)
  number_of_samples <- as.integer(samples_in_cancer %>% filter(cancer == cancer_type) %>% select(n))
  mutational_matrix[, column] <- mutational_matrix[, column] / number_of_samples
}





ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))

#Create a cosine similarity matrix
cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix, cosmic_signatures)

#Fit mut_matrix to cosmic
fit_to_cosmic <- fit_to_signatures(mutational_matrix,cosmic_signatures)

contribution <- fit_to_cosmic$contribution



#Optional filter of contribution 
contribution <- sweep(contribution,2,colSums(contribution),"/")




contribution_DF <- data.frame()

for (i in 1:ncol(contribution)){
  
  col_vector <- sort.default(contribution[, i],decreasing = TRUE)
  top_x_signatures <- head(col_vector,n=6)
  other <- sum(col_vector[!col_vector %in% top_x_signatures])
  append_df <- data.frame("Contribution"=as.vector(top_x_signatures),"Signature"=names(top_x_signatures))
  append_df <- rbind(append_df,data.frame("Contribution"=other,"Signature"=c("other")))
  
  append_df$cancer <- rep(colnames(contribution)[i],nrow(append_df))
  
  contribution_DF <- rbind(contribution_DF,append_df)
  
}


rm(append_df,other,col_vector,top_x_signatures,i,contribution)





#contribution_DF <- plyr::arrange(contribution_DF,Signature)
#contribution_DF <- contribution_DF[order(contribution_DF$Signature, decreasing = T),]



ggplot(data=contribution_DF, aes(x=cancer, y=Contribution,group = Signature)) +
  geom_bar(aes(fill= Signature),stat="identity",color="black",position = "Stack")+
  #geom_col(position = position_stack(reverse = TRUE))+
  geom_text(aes(y=Contribution, label=Signature), vjust=1, 
            color="black", size=2.5,position="Stack")+
  theme_minimal()+
  theme(axis.text.x=element_text(vjust = 0.5,size=10,angle = 90))



  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#Alternative plotting methods
  
  p<-ggplot(data=contribution_DF, aes(x=cancer, y=Contribution,fill=Signature)) +
    geom_bar(stat="identity",color="black")+
    geom_text(aes(label=Signature), vjust=1.6, color="white", size=3.5)+
    theme_bw()+
    theme(axis.text.x=element_text(vjust = 0.5,size=10,angle = 90))
  p
  

  library(plyr)
  # Sort by dose and supp
  df_sorted <- plyr::arrange(contribution_DF,cancer,Signature)
  
  
  df_cumsum <- ddply(df_sorted, "cancer",
                     transform, label_ypos=cumsum(Contribution)-0.5*Contribution)
  
  
  p<-ggplot(data=df_cumsum, aes(x=cancer, y=Contribution,fill=Signature)) +
    geom_bar(stat="identity",color="black")+
    #geom_col(position = position_stack(reverse = TRUE))+
    geom_text(aes(y=label_ypos, label=Signature), vjust=1.6, 
              color="white", size=3.5)+
    theme_minimal()+
    theme(axis.text.x=element_text(vjust = 0.5,size=10,angle = 90))
  p
  
  




#Example pipeline
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

#Load files with DNA_repair proteins
DNA_repair <- read.csv("DNA_repair_extended.csv",sep=";")

#Set reference genome
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"




# MC3 file  ---------------------------------------------------------------------


#Transform MC3 file to VCF files (needed for mutational_patterns)
#only need to be done once.
if (!dir.exists(file.path(paste(main_wd,"/VCF_files",sep="")))){
  
  dir_name <- paste(main_wd,"/VCF_files",sep="")
  
  dir.create(file.path(dir_name))
  MC3_to_VCF()
  
}  
  


#Read MC3 file
if(!file.exists("MC3.rda")){
  
  MC3 <- read_MC3()
  save(MC3,file="MC3.rda")
  
}
load("MC3.rda")


 

#Random samples
#random_samples <- unique(MC3$Tumor_Sample_Barcode)
#random_samples <- sample(random_samples,200)

#high_mutation_samples <- random_samples




#Extract highly mutated samples
#For convenince S and type is how files are named such as pictures or RDA objects
# 1000 <- N=14, 2000 <- N=8 , 3 <- N=14

N <- 14
S <- 3
type <- "SD" #SD or absolute

high_mutation_samples <- get_high_mutations(MC3,type,S)


#Subselction of MC3 

#Select genes associated with DNA repair mechanisms
#MC3 <- MC3 %>% filter(Hugo_Symbol %in% DNA_repair$hgnc_symbol)

#Select samples with high mutation rates
MC3 <- MC3 %>% filter(Tumor_Sample_Barcode %in% high_mutation_samples)



  
# Mutational signatures ---------------------------------------------------

vcf_list <- list.files(path=paste(main_wd,"/VCF_files",sep=""),recursive=TRUE,pattern=".vcf")

#select only highly mutated samples
vcf_list_names <- gsub(".vcf","",gsub("[A-Z0-9]*/","",vcf_list))
vcf_list <- vcf_list[vcf_list_names %in% high_mutation_samples]
vcf_list_names <- gsub(".vcf","",gsub("[A-Z0-9]*/","",vcf_list))

#convert VCF files to GRange objects(takes some time)


if (!dir.exists(file.path(paste(main_wd,"/output",sep="")))){
  
  dir_name <- paste(main_wd,"/output",sep="")
  dir.create(file.path(dir_name))
  setwd(paste(main_wd,"/VCF_files",sep=""))
  all_vcfs <- read_vcfs_as_granges(vcf_list,vcf_list_names,ref_genome)
  
  setwd(dir_name)
  save(all_vcfs,file=paste(S,"_all_vcfs.rda",sep=""))
}  

file_name <- paste(S,"_all_vcfs.rda",sep="")
dir_name <- paste(main_wd,"/output",sep="")
setwd(dir_name)

if(!file.exists(file_name)){
  setwd(paste(main_wd,"/VCF_files",sep=""))
  all_vcfs <- read_vcfs_as_granges(vcf_list,vcf_list_names,ref_genome)
  setwd(dir_name)
  save(all_vcfs,file=paste(S,"_all_vcfs.rda",sep=""))
  
} else {
  setwd(dir_name)
  load(file_name)
  
}

#dir_name <- paste(main_wd,"/output",sep="")
#load(paste(S,"_all_vcfs.rda",sep=""))



#Create a mutational matrix
setwd(dir_name)
file_name <- paste(type,S,"mut_mat.rda",sep="_")
if(!file.exists(file_name)){
  
  mutational_matrix = mut_matrix(all_vcfs,ref_genome)
  save(mutational_matrix,file=file_name)
} else {
  
  load(file_name)
  
}



#Extract cosine similarity
setwd(main_wd)
#Read cosmic signatures
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))


#Just plot number in heatmap.
sampleDF <- data.frame("TCGA_code"=colnames(mutational_matrix))
sampleDF$sample <- c(1:ncol(mutational_matrix))
colnames(mutational_matrix) <- c(1:ncol(mutational_matrix))


#Shorten TCGA-barcode
#colnames(mutational_matrix) <- gsub("TCGA-","",
#                                    gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(mutational_matrix)))

#Create a cosine similarity matrix
cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix, cosmic_signatures)

#Cluster cosmic signatures
hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
#store signatures in new order
cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]

#Plot everything in one pdf.file  PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF PDF
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Pictures/Test")
#pdf(file=paste(S,"_high_mut.pdf"),height=8.27,width =11.69)


#plot cosine_heatmap
print(plot_cosine_heatmap(cos_sim_samples_cosmic, col_order = cosmic_order, cluster_rows = TRUE))



#Cluster samples based on cosine similarity
#N <- 14
sample_cluster <- hclust(dist(cos_sim_samples_cosmic,method="euclidean"),method="complete")
plot_cluster(sample_cluster,N,2.2)




setwd(main_wd)
#plot cluster in heatmap
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



if(FALSE){
  #Manipulate extracted cluster
  #For >1000 cancer. N = 14
  
  
  remove_cluster <- c(5,10,13,14)
  ClusterDF <- ClusterDF %>% filter(!cluster %in% remove_cluster)
  ClusterDF[ClusterDF$cluster == 5 , 2] <- 4  
  ClusterDF[ClusterDF$cluster == 6 , 2] <- 5
  ClusterDF[ClusterDF$cluster == 7 , 2] <- 5
  ClusterDF[ClusterDF$cluster == 8 , 2] <- 6
  ClusterDF[ClusterDF$cluster == 9 , 2] <- 7
  ClusterDF[ClusterDF$cluster == 10 , 2] <- 8
  ClusterDF[ClusterDF$cluster == 12 , 2] <- 9
  ClusterDF[ClusterDF$cluster == 13 , 2] <- 10
  
  plot_cosine_heatmap()
  
  
  
  
}

if(TRUE){
  #Remove small cluster
  keep_cluster <- ClusterDF %>% dplyr::count(cluster)
  keep_cluster <- (keep_cluster[keep_cluster$n > 3 , 1])
  keep_cluster <- keep_cluster$cluster
  
  ClusterDF <- ClusterDF[ClusterDF$cluster %in% keep_cluster ,]
  
  
  
}


rm(TSS2Study,sampleDF,MC3)



#Find signatures in cluster
#Make sure that rownames of cos_sim matches ClusterDF$samples
#get_signature(ClusterDF,cos_sim_samples_cosmic,0.6)


#Main contributing signatures. (Those needed to achieve the threshold 
#cosine simililarity between reconstruction and original)
contributing_signatures <- get_contributing_signatures(ClusterDF,cos_sim_samples_cosmic,
                            mutational_matrix,threshold = 0.90)


library(gridExtra)
library(grid)
library(rowr)

contributing_signatures_DF <- data.frame()

for (i in 1:length(contributing_signatures)){
 
  
  signatures_vector <- paste(contributing_signatures[[i]],collapse=",")
  name_vector <- names(contributing_signatures)[i]
  append_df <- data.frame("Cluster" = name_vector,"Prominent signatures" = signatures_vector)
  
  print(append_df)
  
  contributing_signatures_DF <- rbind(contributing_signatures_DF,append_df)
  
}

colnames(contributing_signatures_DF) <- c("Cluster","Prominent Signatures")


#pdf(file = "contsign.pdf", height = 12, width = 26)
grid.newpage()
grid.table(contributing_signatures_DF,rows = NULL)
#dev.off()


contributing_signatures


###############################################################TRY GETTING RIGHT SIGNATURE 

Cluster_id <- 1
signatures <- contributing_signatures$`cluster 1`

pop_signatures <- function(signatures,Cluster_id){
  
  
  samples_in_cluster <- as.vector(ClusterDF %>% filter(cluster == Cluster_id) %>% select(sample))
  samples_in_cluster <- samples_in_cluster$sample
  Mut_Mat_Cluster <- mutational_matrix[, colnames(mutational_matrix) %in% samples_in_cluster]
  cosmic_in_cluster <- cosmic_signatures[, colnames(cosmic_signatures) %in% signatures]
  

  
  may_i_pop <- function(signature_to_pop,all_signatures){
    
    #With
    signatures_to_use <- cosmic_in_cluster[, colnames(cosmic_in_cluster) %in% all_signatures]
    fit_sign <- fit_to_signatures(Mut_Mat_Cluster,signatures_to_use)
    cos_sim_original_reconstructed <- cos_sim_matrix(Mut_Mat_Cluster, fit_sign$reconstructed)
    cos_sim_original_reconstructed <- as.data.frame(diag(cos_sim_original_reconstructed))
    with_cossim <- cos_sim_original_reconstructed
    
    #Without pop_sign
    signatures_to_use <- all_signatures[!all_signatures == signature_to_pop]
    signatures_to_use <- cosmic_in_cluster[, colnames(cosmic_in_cluster) %in% signatures_to_use]
    fit_sign <- fit_to_signatures(Mut_Mat_Cluster,signatures_to_use)
    cos_sim_original_reconstructed <- cos_sim_matrix(Mut_Mat_Cluster, fit_sign$reconstructed)
    cos_sim_original_reconstructed <- as.data.frame(diag(cos_sim_original_reconstructed))
    without_cossim <- cos_sim_original_reconstructed
    
    
    
    
    p_val <- t.test(with_cossim,without_cossim)$p.value
    print(p_val)
    print(signature_to_pop)
    
    
    
    return(p_val > 0.40)
    
    
  }
  
  
  
  #loop until no signature can be removed without significantly changing cossim
  checked_sign <- 0
  i <- 0
  counter <- 1
  tot_sign <- length(signatures)
  #sucess <- FALSE
  
  while(TRUE){
    
    
    print(paste("#sign",length(signatures),sep="="))
    print(paste("i",i,sep="="))
    counter <- counter + 1
    
    if (may_i_pop(signatures[length(signatures)-i],signatures) == TRUE){
      
      signatures <- signatures[!signatures == signatures[length(signatures)-i]]
      #signatures <- head(signatures,-1)
      checked_sign <- 0
      #tot_sign <- length(signatures)
      
      
    }else{
      i <-  i + 1
      checked_sign <- checked_sign + 1
    }
    
    #print(checked_sign)
    #print(signatures)
    #(length(signatures) == 2 | checked_sign == (tot_sign-1))
    
    if (counter == tot_sign - 2){
      
      break()
    }
    
    
    

    
  }
  
  print(paste("FOUND SIGNATURES:",signatures,sep=""))
  
}

#----------------------------------------------------------------------------------------------------------------




apply(contributing_signatures,function(x)data.frame(paste(x,collapse=",")))

empty_df <- data.frame()

for (item in contributing_signatures){
  

  empty_df <- rbind(empty_df,data.frame("signature"=paste(item,collapse = ",")))
  
  
}

empty_df <- data.frame("Cluster"=names(contributing_signatures),"Signature"=as.vector(empty_df))
empty_df$Cluster <- gsub("cluster","",empty_df$Cluster)


library(formattable)
print(formattable(empty_df,align="l"))




piechart_cancer_cluster(ClusterDF,mfrows=c(round(N/2),2))


rm(list=ls()[! ls() %in% c("ClusterDF","DNA_repair","S","type","ref_genome","main_wd")])

# CNV_SNV plot ------------------------------------------------------------
setwd(main_wd)

#samples <- ClusterDF %>% filter(cluster == 3) %>% select(TCGA_code)
#samples <- as.character(samples$TCGA_code)
#samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples)


#genes <- as.character(DNA_repair$hgnc_symbol)

#################################################################################
#Test MMR cluster 6,7,9
#samples <- ClusterDF %>% filter(cluster %in% c(1,4,8)) %>% select(TCGA_code)
samples <- ClusterDF %>% select(TCGA_code)
#samples <- unique(samples)
samples <- as.character(samples$TCGA_code)
samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples)

genes <- read.table("gen_lista.csv",header=TRUE,sep=";")
#genes <- genes %>% filter(System == "BER")
genes <- as.character(genes$Gen)



plot_CNV_SNV(samples,genes,"High mutation samples",plot_id=FALSE,cluster_names = TRUE)




#Just take random genes 
#samples <- as.character(sample(MC3$Tumor_Sample_Barcode,40))
#samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples)
#genes <- as.character(unique(MC3$Hugo_Symbol[1:40]))
#plot_CNV_SNV(samples,genes,"Random")

####################################################

if(FALSE){
  #Plot and save all cluster
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Pictures/CNV_SNV")
  pdf(file="one_cluster_a_time.pdf",width=16, height=10)
  
  
  for (Cluster in unique(ClusterDF$cluster)){
    
    setwd(main_wd)
    
    samples <- ClusterDF %>% filter(cluster == Cluster) %>% select(TCGA_code)
    samples <- as.character(samples$TCGA_code)
    samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples)
    
    
    genes <- as.character(DNA_repair$hgnc_symbol)
    print(plot_CNV_SNV(samples,genes,paste("Cluster",Cluster,sep=" ")))
    
    
  }
  
  dev.off()
  
  
}


# mRNA Exp ----------------------------------------------------------------
#genes <- read.table("gen_lista.csv",header=TRUE,sep=";")
#genes <- as.character(genes$Gen)
#genes <- unique(genes) #unrefined gene list
mRNA_DF <- mRNA_exp(genes)

plot_mRNA_SNV(samples,genes,cluster_names=TRUE,mRNA_DF)



#Try AID/APOBEC 
#setwd(main_wd)
#samples <- ClusterDF %>% filter(cluster %in% c(1,4,8)) %>% select(TCGA_code)

#Select samples
#samples <- ClusterDF$TCGA_code
#samples <- as.character(samples$TCGA_code)
#samples <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",samples)


#Select genes
#genes <- read.table("gen_lista.csv",header=TRUE,sep=";")
#genes <- genes %>% filter(System == "AID/APO")
#genes <- as.character(genes$Gen)





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

  
#In order to not get grey values
exp[exp >2] <- 2
  


library(ggplot2)
library(reshape2)

exp_melt <- melt(exp)

ggplot(data = exp_melt, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + labs(x="",y="",title="High mutation samples")+
  #scale_fill_gradient(low="green", mid="lightblue", high="red",limits
  scale_fill_gradientn(colours=c("red","lightblue","green"), limits=c(0,2))+
  geom_hline(yintercept=hline_pos)+
  
  theme(
    axis.text.x=element_text(colour="black",angle=90, hjust=1,vjust = 0.5,size=10),
    axis.text.y=element_text(vjust = 0.5,size=6)
  )







dev.off()








# Find if any certain SNVS is common within cluster --------------------------------------
library(data.table)
library(dplyr)
#might just aswell be a seperate script


rm(list=ls()[! ls() %in% c("ClusterDF","main_wd","genes")])
gc()

main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)

library(data.table)


#ClusterDF
#ClusterDF<- read.table("ClusterDF.txt",header=TRUE)


#Need information about AA changes
#MC3_DF <- fread('mc3.v0.2.8.PUBLIC.maf.gz')

#MC3_DF <- MC3_DF %>% select("Hugo_Symbol","Chromosome","Start_Position"
                           # ,"End_Position","Tumor_Sample_Barcode","Variant_Classification"
                           # ,"HGVSp_Short","Strand")


#MC3 <- MC3_DF
load("MC3_ext.rda")



genes <- read.table("gen_lista.csv",header=TRUE,sep=";")
#genes <- genes %>% filter(System == "BER")
genes <- as.character(genes$Gen)

table_to_print <- tibble()

for (i in 1:length(unique(ClusterDF$cluster))){
  

  
  cluster_x <- ClusterDF %>% filter(cluster == i)
  
  cluster_MC3 <- MC3 %>% filter(Tumor_Sample_Barcode %in% cluster_x$TCGA_code)
  common_mutations <- cluster_MC3 %>% dplyr::count(Hugo_Symbol,Chromosome, 
                                                   Start_Position,End_Position,
                                                   Variant_Classification,HGVSp_Short)
  
  common_mutations <- common_mutations[order(common_mutations$n,decreasing = TRUE),]
  test <- head(common_mutations,10)
  test <- test %>% select(n,Hugo_Symbol,HGVSp_Short,Variant_Classification)
  
  
  #Check if any of the commonly mutated genes are in the genes-list
  high_mut_gene  <- intersect(test$Hugo_Symbol,genes)
  

  test <- test %>% filter(Hugo_Symbol %in% high_mut_gene)
  
  print(paste("cluster",i,sep=" "))
  print(high_mut_gene)
  
  test$cluster <- c(rep(i,nrow(test)))
  table_to_print <- rbind(table_to_print,test)
  
  #library(formattable)
  #formattable(test,align="l")
  

}

library(formattable)
library(kableExtra)
colnames(table_to_print) <-  c("Number of occurrences","Hugo Symbol","Amino acid change","Variant Classification","cluster")               

print(formattable(table_to_print,align="l"))
#print(kable_styling(kable(table_to_print)))


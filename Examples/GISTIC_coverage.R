#GISTIC coverage
#Percentage of genes mutated in each cluster
#CNV_coverage
rm(list=ls())
gc()
#Library
library(dplyr)
library(ggplot2)
library(data.table)

#Set wd to folder containing all data
main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/GISTIC")
CNV_files <- list.files(pattern = "all_thresholded.by_genes",recursive = TRUE)

#Open and read all treshold files
#thresholded_files <- list.files(pattern = "all_thresholded.by_genes",recursive = TRUE)

#treshold_files_df <- lapply(treshold_files, function(x) {fread(file = x, header = T,stringsAsFactors = FALSE)})

cancer_file_name <- CNV_files[[1]]


coverage_list <- list()
coverage_DF <- data.frame()

for (cancer_file_name in CNV_files){
  
  
  cancer_type <- gsub("..*org_","",gsub("-T..*$","",cancer_file_name))
  print(paste("Starting on",cancer_type,sep=" "))
  cancer_file <- fread(file = cancer_file_name, header = T,stringsAsFactors = FALSE)
  
  
  
  cancer_file <- cancer_file[, 4:ncol(cancer_file)]
  
  count_CNV <- apply(cancer_file, 2, function(x) {as.data.frame(count(x))})
  

  
  get_percentage <-function(CNV){
    
    total_sum <- sum(CNV$freq)
    tot_CNV   <- CNV %>% filter(x != 0)
    tot_CNV   <- sum(tot_CNV$freq)
    percent_with_CNV <- tot_CNV / total_sum
    
    return(percent_with_CNV)
    
  }
  
  
  
  percentage_CNV <- sapply(count_CNV,get_percentage)
  
  coverage_list[[cancer_type]] <- percentage_CNV
  
  append_DF <- data.frame("Percentage"=as.numeric(percentage_CNV),"Sample"=names(percentage_CNV),
                          "Cancer" = rep(cancer_type,length(percentage_CNV)))

  coverage_DF <- rbind(coverage_DF,append_DF)
  
  
}













box_p <- ggplot(coverage_DF, aes(x=reorder(Cancer,Percentage,FUN = median), y=Percentage)) +
  geom_boxplot() +
  stat_boxplot(geom ='errorbar')+
  #theme_bw()+
  labs(title="",x="Cancer", y = "Part of genome with CNV:s")+
  #geom_boxplot(coef = 6) +
  geom_boxplot(outlier.colour=NA,outlier.size = 0,outlier.alpha = 0) +
  #coord_cartesian(ylim = c(0, 3500))+
  theme_minimal()+
  theme(axis.text.x=element_text(vjust = 0.5,size=rel(1.5),angle = 90),
        axis.text.y=element_text(size=rel(1.5)),
        axis.title=element_text(size=rel(1.5)))


box_p

#Remove points
library(plotly)
box_p <- plotly_build(box_p)

for(i in 1:length(box_p$x$data)) {
  box_p$x$data[[i]]$marker$opacity = 0
  
}

box_p












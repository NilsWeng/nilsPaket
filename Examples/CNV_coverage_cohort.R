#CNV_coverage
rm(list=ls())
gc()
#Library
library(dplyr)
library(ggplot2)


#Set wd to folder containing all data
main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)

setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/GISTIC")
CNV_files <- list.files(pattern = "focal_input.seg",recursive = TRUE)

#Open and read all treshold files
treshold_files <- list.files(pattern = "sample_cutoffs",recursive = TRUE)
treshold_files_df <- lapply(treshold_files, function(x) {read.table(file = x, header = T)})
threshold_df <- do.call("rbind", lapply(treshold_files_df, as.data.frame)) 
rm(treshold_files,treshold_files_df)


#Tot length Excluding Y-chromosome (Chr1-X) Ask malin
tot_hg19_length <- 3036303846






#Function for finding coverage of genome
CN_coverage <- function(sample){
  
  #Deep amp/del filtering
  #sample_id <- as.character(unique(sample$Sample))
  #high_thresh <- as.numeric(threshold_df %>% filter(Sample == sample_id) %>% select(High))
  #low_thresh <- as.numeric(threshold_df %>% filter(Sample == sample_id) %>% select(Low))
  #CNV_regions <- sample %>% filter(sample$Seg.CN >= high_thresh  | sample$Seg.CN <= low_thresh  )
  
  
  CNV_regions <- sample %>% filter(abs(sample$Seg.CN) >= 0.1)
  seg_length <- sum(CNV_regions$End.bp-CNV_regions$Start.bp)
  CNV_coverage <- seg_length / tot_hg19_length
  
  return(CNV_coverage)
  
  
}






#Loop over each cancer_type and count length of CN changed regions

coverage_DF <- data.frame()

for (cancer_CNV_file in CNV_files) {
  
  cancer <- gsub("gdac.broadinstitute.org_","",gsub("-..*$","",cancer_CNV_file))
  CNV_cancer <- read.table(cancer_CNV_file,header=TRUE)
  CNV_samples <- split.data.frame(CNV_cancer,CNV_cancer$Sample)
  
  
  
  
  Cancer_CNV_coverage <- as.data.frame(lapply(CNV_samples,CN_coverage))
  append_df <- data.frame("sample"=names(Cancer_CNV_coverage),"CNV_coverage"=as.numeric(Cancer_CNV_coverage),
                          "cancer" = rep(cancer,length(names(Cancer_CNV_coverage))))
  
  
  coverage_DF <- rbind(coverage_DF,append_df)
  
  
}


avarage_coverage <- coverage_DF %>%
  group_by(cancer) %>%
  dplyr::summarize(Mean = mean(CNV_coverage, na.rm=TRUE))


#Boxplot of coverage for each cancer type


box_p <- ggplot(coverage_DF, aes(x=reorder(cancer,CNV_coverage,FUN = median), y=CNV_coverage)) +
  geom_boxplot() +
  stat_boxplot(geom ='errorbar')+
  #theme_bw()+
  labs(title="",x="Cancer", y = "Percentage of genome with abs(CN) > 0.1")+
  #geom_boxplot(coef = 6) +
  geom_boxplot(outlier.colour=NA,outlier.size = 0,outlier.alpha = 0) +
  #coord_cartesian(ylim = c(0, 3500))+
  theme_minimal()+
  theme(axis.text.x=element_text(vjust = 0.5,size=8,angle = 90))


box_p

#Remove points
library(plotly)
box_p <- plotly_build(box_p)

for(i in 1:length(box_p$x$data)) {
  box_p$x$data[[i]]$marker$opacity = 0
  
}

box_p











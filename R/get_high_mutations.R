#' Extract high mutation samples
#'
#' Extract high mutation samples
#'
#' @param MC3_files MC3 data frame created by read_MC3
#' @param type sets treshold type , either "absolute" or "SD"
#' @param treshold numerical treshold value either absolute number or number of SD:s 
#'
#' @return None
#' @export
#' @examples
#' read_GISTIC(MC3,type="SD",3)
#' @import dplyr
#' @import plyr



get_high_mutations <- function(MC3_file,type,treshold){
  
  
  counts <- as.data.frame(table(MC3_file$Tumor_Sample_Barcode))
  counts$TSS.Code <- gsub("TCGA-","",gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",counts$Var1))
  
  #Get cancertype for each sample
  TSS2Study_DF <- read.table("TSS2Studyabb.txt",header=TRUE)
  counts <- join(counts,TSS2Study_DF, by="TSS.Code")
  
  
  high_mut_samples_vector <- c()
  
  for (cancer in as.vector(unique(counts$Study.Abbreviation))){
    
    
    counts_cancer <- counts %>% filter(Study.Abbreviation == cancer)
    
    
    if(type == "absolute"){
      
      
      
      high_mut_samples <- as.vector(unlist(counts_cancer %>% filter(Freq > treshold) %>% select(Var1)))
      
    }
    
    
    if(type == "SD"){
      
      
      
      median <- as.numeric(median(counts_cancer$Freq))
      S <- as.numeric(sd(counts_cancer$Freq))
      treshold1 <- as.numeric(median + treshold*S)
      
      high_mut_samples <- as.vector(unlist((counts_cancer %>% filter(Freq > treshold1) %>% select(Var1))))
      
      
    }
    
    high_mut_samples_vector <- c(high_mut_samples_vector,high_mut_samples)
    
    
  }
  
  
  
  high_mut_samples_vector
  
  return(high_mut_samples_vector)  
  
}




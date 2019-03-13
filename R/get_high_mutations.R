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




get_high_mutations <- function(MC3_file,type,treshold){
  
  
  counts <- MC3_file %>% count(Tumor_Sample_Barcode)
  
  
  if(type == "absolute"){
    
    
    
    high_mut_samples <- counts %>% filter(n > treshold) %>% select(Tumor_Sample_Barcode)
    
  }
  
  
  if(type == "SD"){
    
    
    
    median <- median(counts$n)
    S <- sd(counts$n)
    treshold <- median + treshold*S
    
    high_mut_samples <- counts %>% filter(n > treshold) %>% select(Tumor_Sample_Barcode)
    
    
    
  }
  
  high_mut_samples <- as.vector(high_mut_samples$Tumor_Sample_Barcode)
  
  return(high_mut_samples)  
  
}




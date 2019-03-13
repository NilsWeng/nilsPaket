#' Extract high mutation samples
#'
#' Extract high mutation samples
#'
#' @param treshold numerical treshold value setting number of standard deviations
#'
#' @return None
#' @export
#' @examples
#' high_mut(3)

#' @import dplyr



high_mut <- function(treshold){
  
  
  counts <- MC3 %>% count(Tumor_Sample_Barcode)
  median <- as.numeric(median(counts$n))
  S <- as.numeric(sd(counts$n))
  treshold <- as.numeric(median + treshold*S)
  
  high_mut_samples <- counts %>% filter(n > treshold) %>% select(Tumor_Sample_Barcode)
  
  return(as.vector(high_mut_samples$Tumor_Sample_Barcode))
  
}
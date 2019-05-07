#' Extract high mutation samples
#'
#' Extract high mutation samples
#'
#' @param treshold What top percentile should be taken from each cancer type
#'
#' @return high_cnv
#' @export
#' @examples
#' read_GISTIC(MC3,type="SD",3)
#' @import dplyr
#' @import plyr
#' @import data.table




get_high_CNV_samples <- function(threshold){
  
  gistic_files <- list.files(pattern="all_thresholded.by_genes",recursive = TRUE)
  gistic_samples <- c()
  high_CNV <- c()
  
  for (file_name in gistic_files){
    
    gistic_file <- fread(file_name,header=TRUE,stringsAsFactors = FALSE)
    gistic_file <- gistic_file[, 4:ncol(gistic_file)]
    
    
    
    count_CNV <- apply(gistic_file, 2, function(x) {as.data.frame(table(x))})
    
    
    
    get_percentage <-function(CNV){
      
      total_sum <- sum(CNV$Freq)
      tot_CNV   <- CNV %>% filter(x != 0)
      tot_CNV   <- sum(tot_CNV$Freq)
      percent_with_CNV <- tot_CNV / total_sum
      
      return(percent_with_CNV)
      
    }
    
    
    
    percentage_CNV <- sapply(count_CNV,get_percentage)
    
    #mean_CNV <- mean(percentage_CNV)
    #S_CNV <- sd(percentage_CNV)
    
    
    percentage_CNV <- sort(percentage_CNV,decreasing = FALSE)
    
    
    high_CNV <- c(high_CNV,names(which(percentage_CNV > quantile(percentage_CNV,threshold))))
    
    
    
  }
  
  return(high_CNV)
  
  
}
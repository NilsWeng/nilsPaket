#' Read GISTIC file
#'
#' Read Gistic file
#'
#' @param files2read Vector of files to read
#' @param selectProteins Vector of proteins names to only select from Defaults to TRUE
#'
#' @return None
#' @export
#' @examples
#' read_GISTIC(files,DNArepair_proteins)

#' @import dplyr
#' @import data.table



read_GISTIC <- function(files2read,selectProteins){
  
  
  read_file <- function(files){
    
    
    
    Table <- fread(files,header=TRUE,stringsAsFactors = FALSE)
    #Filter only for genes involved in DNA-repair
    Table <- Table %>% filter(Table$`Gene Symbol` %in% DNA_repair$hgnc_symbol)
    return(Table)
    
    
  }
  
  
  
  dataset <- do.call("cbind",lapply(files2read ,read_file))
  dataset <- dataset[, !duplicated(colnames(dataset))]
  
  return(dataset)
  



}





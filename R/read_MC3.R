#' Read MC3 file
#'
#' Read MC3 file
#'
#' @param selectProteins Vector of genenames to make a subselction from
#'
#' @return MC3_DF
#' @export
#' @examples
#' read_MC3(DNArepair_genes)

#' @import dplyr
#' @import data.table


read_MC3 <- function(genes2select){
  
  MC3_DF <- fread('mc3.v0.2.8.PUBLIC.maf.gz')
  MC3_DF <- MC3_DF %>% select("Hugo_Symbol","Chromosome","Start_Position"
                              ,"End_Position","Tumor_Sample_Barcode","Variant_Classification")
  
  if(missing(genes2select)){
    
    MC3_DF <- MC3_DF
    
  } else {
    
    MC3_DF <- MC3_DF[MC3_DF$Hugo_Symbol %in% genes2select ,]
    
  }
  
  return(MC3_DF)
  
}
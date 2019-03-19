#' mRNA_exp
#'
#' Given some genes get mRNA expression for all samples in TCGA cohort
#'
#' @param genes Genes to make a subselection from (to large otherwise)
#'
#' @return mRNA_DF
#' @export
#' @examples
#' mRNA_exp(c("ACC","XX","YY"))

#' @import dplyr
#' @import data.table

mRNA_exp <- function(genes){
  
  setwd(main_wd)
  file_list <- list.files(pattern="_normalized__data.data.txt",recursive = TRUE)
  print("reading mRNA files")
  
  read_mRNASeq <- function(files){
    
    
    
    Table1 <- fread(files,header=TRUE,stringsAsFactors = FALSE)
    #Filter only for genes in function argument.
    Table1$`Hybridization REF` <- gsub("\\|.*$","",Table1$`Hybridization REF`)
    Table1 <- Table1 %>% filter(Table1$`Hybridization REF` %in% genes)
    
    
    return(Table1)
    
    
  }
  
  #file_list <- file_list[1:3]
  
  mRNA_DF <- do.call("cbind",lapply(file_list ,read_mRNASeq))
  mRNA_DF <- mRNA_DF[, !duplicated(colnames(mRNA_DF))]
  rownames(mRNA_DF) <- as.vector(mRNA_DF$`Hybridization REF`)
  mRNA_DF<- mRNA_DF[, !colnames(mRNA_DF)=="Hybridization REF"]
  colnames(mRNA_DF) <-  gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(mRNA_DF))
  
  
  setwd(main_wd)
  #write.table(dataset,"PANCAN_mRNA_expression.txt")
  #save(dataset,file="DNArep_mRNA_exp.rda")
  return(mRNA_DF)
}
#' commonly_mut_genes
#' 
#' Find genes that are commonly mutated from ClusterDF
#' 
#' 
#' @param ClusterDF  ClusterDF
#' @export
#' @return Table with which genes in which cluster
#' @import kableExtra
#' @import dplyr
#'
#' @examples commonly_mut_genes(ClusterDF)


commonly_mut_genes <- function(ClusterDF){
  
  setwd(main_wd)
  load("MC3_ext.rda")
  genes <- read.table("gen_lista.csv",header=TRUE,sep=";")
  genes <- as.character(genes$Gen)
  
  
  
  table_to_print <- tibble()
  
  for (i in 1:length(unique(ClusterDF$cluster))){
    
    
    
    cluster_x <- ClusterDF %>% filter(cluster == i)
    
    cluster_MC3 <- MC3 %>% filter(Tumor_Sample_Barcode %in% cluster_x$TCGA_code)
    common_mutations <- cluster_MC3 %>% dplyr::count(Hugo_Symbol,Chromosome, 
                                                     Start_Position,End_Position,
                                                     Variant_Classification,HGVSp_Short)
    
    common_mutations <- common_mutations[order(common_mutations$n,decreasing = TRUE),]
    test <- head(common_mutations,30)
    test <- test %>% select(n,Hugo_Symbol,HGVSp_Short,Variant_Classification)
    
    
    #Check if any of the commonly mutated genes are in the genes-list
    high_mut_gene  <- intersect(test$Hugo_Symbol,genes)
    
    
    test <- test %>% filter(Hugo_Symbol %in% high_mut_gene)
   
    
    test$cluster <- c(rep(i,nrow(test)))
    table_to_print <- rbind(table_to_print,test)

    
  }
  
  
  library(kableExtra)
  print(formattable(table_to_print,align="l"))
  
  return(table_to_print)
}
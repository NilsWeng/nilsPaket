#' get_stat_as_matrix 
#'
#' Function that given vector of samples and genes and a mRNA expression table(from mRNA_exp function) returns a matrix with exp/avg exp
#' @param samples vector of samples
#' @param genes vector of genes
#' @param DF DF from mRNA_exp function
#' @export
#' @examples
#' get_stat_as_matrix(samples,genes,DF)

#' @import dplyr 
#' @return mRNA_exp.m



get_stat_as_matrix <- function(samples,genes,DF){
  #Make sure to only use samples that have mRNA exp (DF)
  genes <- genes[genes %in% rownames(DF)]
  samples <- samples [samples %in% colnames(DF)]
  
  genes <- unique(genes)
  samples <- unique(samples)
  
  mRNA_exp.m <-  matrix(nrow=length(samples),ncol = length(genes))
  colnames(mRNA_exp.m) <- genes
  rownames(mRNA_exp.m) <- samples
  DF <- DF[match(genes,rownames(DF)),]
  
  
  cancer_vector <- get_cancer_type(colnames(DF))
  cancerDF <- data.frame("Cancer"=cancer_vector,"sample"=colnames(DF))
  
  for (sample in samples){
    
    exp <- as.numeric(DF[, colnames(DF) == sample])
    
    cancer_type <- as.character(cancerDF[cancerDF$sample == sample ,1])
    cancer_samples <- cancerDF %>% filter(Cancer==cancer_type) %>% select(sample)
    cancer_samples <- as.character(cancer_samples$sample)
    cancer_samples <- cancer_samples[!cancer_samples == sample]
    cancer_exp <- DF[, colnames(DF) %in% cancer_samples]
    cancer_exp <- as.matrix(cancer_exp)
    mode(cancer_exp) <- "numeric"
    
    stat_vector <- c()
    type_vector <- c()
    
    for (row in 1:nrow(cancer_exp)){
      stat <- wilcox.test(as.numeric(cancer_exp[row ,]),as.numeric(exp[row]))
      stat_vector <- c(stat_vector,stat$p.value)
      
      if (stat$p.value  < 0.2){
        
        diff <- as.numeric(exp[row]) - mean(as.numeric(cancer_exp[row ,]))
        
        if (diff > 0){
          type_vector <- c(type_vector,"Increase")
          
        }else{
          type_vector <- c(type_vector,"Decrease")
        }
        
        
        
        
      }else{
        
        type_vector <- c(type_vector,"Not significant")
        #WStuck here we got type_vector and p_vector for one sample.
        
      }
      
      
    }
  
    
 
    
  }
  
  
  
  return(mRNA_exp.m)
  
  
}

#' get_contributing_signatures
#'
#' Function that finds the minimun number of signatures needed recreated mutational profiles of each cluster to a certain threshold
#' However cant deal with clusters of size = 1. These are removed.
#'
#' @param ClusterDF Dataframe with cluster information 
#' @param cos_sim_matrix cosine simililarity matrix 
#' @param treshold Treshold for needed cosine similarity between reconstructed and original mutational matrix
#' @param mutational_matrix Mutational matrix from mut_mat()
#'
#' @return DF with main contributing signatures in each cluster
#' @export
#' @examples
#' get_signatre(Cluster,cos_sim_matrix,mutational_matrix,treshold=0.95)

#' @import dplyr









get_contributing_signatures_significance <- function(ClusterDF,cos_sim_samples_cosmic,mutational_matrix,threshold) {
  
  
  
  contributing_signatures <- function(cluster_id){
    
    
    
    clust_samples <- ClusterDF %>% filter(cluster == cluster_id) %>% select(sample)
    
    
    
    clust_samples <- clust_samples$sample
    clust_cossim <- cos_sim_samples_cosmic[rownames(cos_sim_samples_cosmic) %in% clust_samples ,]
    top_sign <- names(sort(colMeans(clust_cossim),decreasing = TRUE))
    
    
    #Reconstruc mutational profile using top signature
    clust_mut_mat <- mutational_matrix[, colnames(mutational_matrix) %in% clust_samples]
    
    #Function needed for loop
    get_cosim <- function(counter){
      
      signatures_to_use <- top_sign[1:counter]
      signatures_to_use <- cosmic_signatures[, colnames(cosmic_signatures) %in% signatures_to_use]
      fit_sign <- fit_to_signatures(clust_mut_mat,signatures_to_use)
      
      #How well are the mutational profiles reconstructed using the top i signatures? 
      cos_sim_original_reconstructed <- cos_sim_matrix(clust_mut_mat, fit_sign$reconstructed)
      cos_sim_original_reconstructed <- as.data.frame(diag(cos_sim_original_reconstructed))
      return(cos_sim_original_reconstructed)
      
    }
    
    
    
    
    cosine_similiarity_DF <- data.frame()
    
    for (i in 2:length(top_sign)){
      
      
      
      #get_cosim(i)
      #get_cosim(i+1)
      
      p_val <- t.test(get_cosim(i),get_cosim(i+1))$p.value
      

      
      if(p_val >= 0.05){
        
        prominent_signatures <- top_sign[1:i]
        return(prominent_signatures)
        break
        
      }
      
      
      
    }

    
    
    
  }
  
  
  
  
  #Remove clusters with only one sample
  cluster_ids <- names(table(ClusterDF$cluster)[(table(ClusterDF$cluster) > 1)])
  cluster_ids <- as.vector(cluster_ids)
  
  
  prominent_signatures <- lapply(cluster_ids,contributing_signatures)
  names(prominent_signatures) <- paste("cluster ",cluster_ids,sep="")
  prominent_signatures
  
  
  
  
  
}

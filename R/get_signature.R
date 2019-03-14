#' get_signature
#'
#' Function that return what clusters contain what signatures(takes mean)
#'
#' @param Cluster 
#' @param cos_sim_matrix 
#' @param Treshold Treshold for colmean over cos_sim_similarity
#'
#' @return DF with signatures and their colmean of cosSim
#' @export
#' @examples
#' get_signatre(Cluster,cos_sim_matrix,Treshold=0.7)

#' @import dplyr


get_signature <- function(ClusterDF,cos_sim_matrix,treshold){
  
  
  signatures_in_cluster <- function(cluster_id){
    
    number_samples <- length(which(ClusterDF$cluster == cluster_id))
    samples <- ClusterDF %>% filter(cluster == cluster_id) %>% select(sample)
    samples_cosine <- cos_sim_matrix[rownames(cos_sim_matrix) %in% samples$sample ,]
    
    if(number_samples > 1){
      
      samples_cosine_mean <- colMeans(samples_cosine)
      
      sig_signatures <- samples_cosine_mean[samples_cosine_mean > treshold]
      sig_signatures <- sort(sig_signatures,decreasing = TRUE)
      
    } else {
      #To handle clusters with only one sample
      samples_cosine_mean <- samples_cosine
    }
    
    sig_signatures <- samples_cosine_mean[samples_cosine_mean > treshold]
    sig_signatures <- sort(sig_signatures,decreasing = TRUE)
    
    return(sig_signatures)
    
  }
  
  
  signatures <- lapply(unique(ClusterDF$cluster),signatures_in_cluster)
  names(signatures) <- unique(ClusterDF$cluster)
  signatures
  
  
  
  
}

#' plot_cluster_in_cosine
#' 
#' Plot cosine similiarity with clustered samples
#' 
#' 
#' @param cluster Clustered samples from hclust
#' @param cos_sim_matrix cosine_similarity matrix from cos_sim_matrix(mutational patterns)
#' @param N Number of clusters to create from cutree.
#'
#' @export
#' @return DF with cluster identity for each sample
#' @import ggplot2
#' @import grid
#' @import zoo
#' @importFrom reshape2 melt
#'
#' @examples plot_cluster_in_cosine(cluster,cos_sim_matrix,8)


plot_cluster_in_cosine <- function(cluster,cos_sim_matrix,N){
  
  cluster_groups <- cutree(cluster,k=N)
  sample_order <- as.vector(unique(cluster_groups[cluster$order]))
  cluster_order <- as.vector(unique(cluster_groups[cluster$order]))
  cluster_groups <- data.frame("sample"=names(cluster_groups),"cluster"=as.vector(cluster_groups))
  cluster_groups <- cluster_groups[order(match(cluster_groups$cluster, sample_order)), ]
  
  #cluster_groups <- cluster_groups[order(cluster_groups$cluster) ,]
  rownames(cluster_groups) <- NULL
  
  sample_order = as.vector(cluster_groups$sample)
  
  hline_pos <- c()
  for (cluster1 in unique(cluster_groups$cluster)){
    
    pos <- tail(which(cluster_groups$cluster == cluster1),n=1)
    pos <- pos + 0.5
    
    hline_pos <- c(hline_pos,pos)
    
  }
  
  
  Cosine.sim = NULL
  Signature = NULL
  Sample = NULL
  x = NULL
  y = NULL
  xend = NULL
  yend = NULL
  # melt
  cos_sim_matrix.m = melt(cos_sim_matrix)
  # assign variable names
  colnames(cos_sim_matrix.m) = c("Sample", "Signature", "Cosine.sim")
  
  

  
  
  # change factor levels to the correct order for plotting
  cos_sim_matrix.m$Signature = factor(cos_sim_matrix.m$Signature, levels = cosmic_order)
  cos_sim_matrix.m$Sample = factor(cos_sim_matrix.m$Sample, levels = sample_order)
  # plot heatmap
  heatmap = ggplot(cos_sim_matrix.m, aes(x=Signature, y=Sample, fill=Cosine.sim, order=Sample)) + 
    geom_tile(color = "white") +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Cosine \nsimilarity", limits = c(0,1)) +
    geom_hline(yintercept=hline_pos)+
    coord_cartesian(clip = 'off') +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=12),
          #axis.text.y = element_text(size=6),
          axis.title.y=element_blank(), # These 3 followin are new
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          
          plot.margin=unit(c(0.5,0.5,0.5,3),"cm")) +
    labs(x=NULL, y=NULL)
  
  print(heatmap)
  
  
  
 
  
  library(grid)
  
  top <- 0.98
  bot <- 0.085
  
  units <- (top-bot)/length(sample_order)
  hline_pos <- c(0,hline_pos)
  hline_pos <- rollmean(hline_pos,k=2)
  hline_pos <- units * hline_pos
  hline_pos <- hline_pos + 0.085
  
  
  my_text <- paste(rep("Cluster",N),cluster_order,sep=" ")
  
  grid.text(my_text,x=0.05, y=hline_pos[1:14],gp=gpar(col="firebrick", fontsize=10, fontface="bold"))
  
  

  return(cluster_groups)
  
}

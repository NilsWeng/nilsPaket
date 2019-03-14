#' plot_cluster 
#'
#' plot_cluster in order to visualise cut-tree
#'
#' @param cluster  Cluster to cut and visualise
#' @param N Number of clusters to be created
#' @param y_pos Define ypos where cluster identity is plotted (2.2)
#' @export
#' @examples
#' plot_cluster(cluster,7,2.2)

#' @importFrom zoo rollmean


plot_cluster <- function(cluster,N,y_pos){
  
  cluster_groups <- cutree(sample_cluster,k=N)
  print(plot(sample_cluster, cex = 0.5,labels=FALSE,xlab="",sub="",ylab="",axes=F,main=""))
  print(rect.hclust(sample_cluster, k = N, border = 2:5))
  
  rect.hclust.labels <- function(cluster,k){
    
    X <- table(cluster_groups)[unique(cluster_groups[cluster$order])]
    m <- c(0, cumsum(X))
    xpos <- rollmean(m,k=2)
    ypos <- rep(y_pos,length(X))
    text(xpos, ypos, names(X), adj=c(0,0), pos=3, col="black", cex=0.75)
    
  }
  
  print(rect.hclust.labels(sample_cluster,k=N))
  
  
}

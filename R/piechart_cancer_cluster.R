#' piechart_cancer_cluster
#'
#' Plot piechart over cluster cancer dist
#' @param ClusterDF blahblah
#' @param mfrows Subplot layout(depends on cluster)
#' @export
#' @examples
#' piechart_cancer_cluster(DF,c(2,3))

#' @import ggplot2
#' @import dplyr





piechart_cancer_cluster <- function(ClusterDF,mfrows=c(3,2)){
  

  setwd(main_wd)
  color_df <- read.table("color_df.txt",header=TRUE)
  
  
  par(mfrow=mfrows)
  par(mar = c(1,0,1,0))
  
  
  for (cluster1 in unique(ClusterDF$cluster)){
    
    
    clust1 <- ClusterDF %>% filter(cluster == cluster1)
    cancer_DF <- clust1 %>% dplyr::count(cluster, Study.Abbreviation)
    cancer_DF <- left_join(cancer_DF,color_df,by="Study.Abbreviation")
    
    
    
    pie(cancer_DF$n,labels = cancer_DF$Study.Abbreviation ,
        col = cancer_DF$colour,main=paste("cluster ",cluster1,"  ( n=", nrow(clust1) ,")",sep=""))
    
    
    
  }
  
  
  par(mfrow=c(1,1))
  legend("right", legend=color_df$Study.Abbreviation, cex=0.8,fill=color_df$colour,xpd="NA")
  
  
}
  




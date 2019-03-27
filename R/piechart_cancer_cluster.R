#' piechart_cancer_cluster
#'
#' Plot piechart over cluster cancer dist
#' @param ClusterDF blahblah
#' @export
#' @examples
#' piechart_cancer_cluster(DF)

#' @import ggplot2
#' @import dplyr



piechart_cancer_cluster <- function(ClusterDF){
  

  
  
  cancer_DF <- ClusterDF %>% dplyr::count(cluster, Study.Abbreviation)
  color_vec <- rainbow(length(unique(ClusterDF$Study.Abbreviation)))
  name_vec  <- as.vector(unique(ClusterDF$Study.Abbreviation))
  color_DF  <- data.frame("Cancer"=name_vec,"color"=color_vec)
  
  #x1 <- as.data.frame(cancer_DF %>% filter(cluster==3))
  #x <- as.data.frame(cancer_DF %>% filter(cluster==4))

  
  
  
  par(mfrow=c(3,2))
  par(mar = c(1,0,1,0))

  for (cluster1 in unique(ClusterDF$cluster)){
    
    x <- as.data.frame(cancer_DF %>% filter(cluster==cluster1))
    x <- x[order(x$Study.Abbreviation),]
    
    
    cancer <- left_join(data.frame("Cancer"=x$Study.Abbreviation),color_DF,by="Cancer")
    color <- as.vector(cancer$color)
    #color <- color_DF %>% filter(Cancer %in% x$Study.Abbreviation) 
    
    #color <- as.vector(color[order(color$Cancer),2])
    
    pie(x$n,labels = x$n, col=color,main=paste("Cluster",cluster1,sep=" "))
    #legend("topright", legend=x$Study.Abbreviation, cex=0.5,fill=color)
    
   
    
    
  }
  
  #plot.new()
  #par(mfrow=c(1,1))
  #legend(0,0, legend = c("0", "1","2","3"),lty = c(2,3,4,5),title = "Cluster",xpd="NA")
  
  #plot.new()
  #par(xpd=NA)
  #legend(locator(1), legend=as.numeric(levels(factor(mtcars$cyl))), pch=19, col= as.numeric(levels(factor(mtcars$cyl))) )
  #legend("topleft", legend=color_DF$Cancer, pch=19,cex=1, col= color_DF$color)
  
  
  par(mfrow=c(1,1))
  legend("right", legend=color_DF$Cancer, cex=0.8,fill=color_DF$color,xpd="NA")
  #plot.new()
  #par(xpd=TRUE)
  #legend("center",legend = color_DF$Cancer,fill=color_DF$color, cex=1)
  #par(xpd=FALSE)
  
  
  
}
  




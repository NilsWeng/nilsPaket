
rm(list=ls())
gc()

main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data")
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
cosmic_signatures <- as.matrix(read.table("cosmic_signatures_extended.txt",header=TRUE))

library(reshape2)


setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Data/mut_matrix")

load("GRange_cohort.rda")


setwd("C:/Users/Nils_/OneDrive/Skrivbord/Main/Pictures/Test")


for (i in 5:10){
  
  
  
  pdf(file=paste(names(GRange_vcf)[i],"_.pdf"),height=8.27,width =11.69)
  
  
  test <- GRange_vcf[[i]]
  l <- split(test,names(test))
  
  print("working on mut mat")
  mutational_matrix <- mut_matrix(l,ref_genome)
  print("done")
  
  plot_mat <- mutational_matrix
  plot_mat <- sweep(plot_mat,2,colSums(plot_mat),"/")
  
  sample_order <- hclust(dist(t(plot_mat)),method="complete")$order
  sample_order <- colnames(plot_mat)[sample_order]
  
  plot_mat <- plot_mat[, match(sample_order ,colnames(plot_mat))]
  
  melt_plot_mat <- melt(plot_mat)
  
  
  p<-ggplot(data=melt_plot_mat,aes(x=Var2, y=Var1, fill=value))+
    geom_tile(color = "gray")+
    scale_fill_gradient(low = "white", high = "blue",limits = c(0,0.2))
  
  
  print(p)
  
  
  
  #Split into high and low mutation samples
  select_high <- which(colSums(mutational_matrix) > 2 * median(colSums(mutational_matrix)))
  select_low  <- which(colSums(mutational_matrix) < 2 * median(colSums(mutational_matrix)))
  
  mutational_matrix_high <- mutational_matrix[, select_high]
  mutational_matrix_low <- mutational_matrix[, select_low]
  
  
  #Cluster rows based on mutational_matrix
  sample_order <- hclust(dist(t(mutational_matrix_high)),method="complete")$order
  sample_order <- colnames(mutational_matrix_high)[sample_order]
  mutational_matrix_high <- mutational_matrix_high[, match(sample_order ,colnames(mutational_matrix_high))]
  
  
  
  #Create a cosine similarity matrix
  cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix_high, cosmic_signatures)
  
  #Fit mut_matrix to cosmic
  fit_to_cosmic <- fit_to_signatures(mutational_matrix_high,cosmic_signatures)
  
  #Filter signatures to plot
  select <- which(rowSums(fit_to_cosmic$contribution)/sum(rowSums(fit_to_cosmic$contribution)) > 0.01)
  
  #Plot
  print(plot_contribution(fit_to_cosmic$contribution[select ,],cosmic_signatures[,select],coord_flip= TRUE,mode="absolute"))
  
  
  #Same for low
  
  
  #Cluster rows based on mutational_matrix
  sample_order <- hclust(dist(t(mutational_matrix_low)),method="complete")$order
  sample_order <- colnames(mutational_matrix_low)[sample_order]
  mutational_matrix_low <- mutational_matrix_low[, match(sample_order ,colnames(mutational_matrix_low))]
  
  
  
  #Create a cosine similarity matrix
  cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix_low, cosmic_signatures)
  
  #Fit mut_matrix to cosmic
  fit_to_cosmic <- fit_to_signatures(mutational_matrix_low,cosmic_signatures)
  
  #Filter signatures to plot Hur smart är det ens att göra såhär? Missar ju sällsynta men kanske betydande 
  select <- which(rowSums(fit_to_cosmic$contribution)/sum(rowSums(fit_to_cosmic$contribution)) > 0.00)
  
  #Plot
  
  print(plot_contribution(fit_to_cosmic$contribution[select ,],cosmic_signatures[,select],coord_flip= TRUE,mode="absolute"))
  
  dev.off()

  
  
  
}



plot_mat <- mutational_matrix
plot_mat <- sweep(plot_mat,2,colSums(plot_mat),"/")

sample_order <- hclust(dist(t(plot_mat)),method="complete")$order
sample_order <- colnames(plot_mat)[sample_order]

plot_mat <- plot_mat[, match(sample_order ,colnames(plot_mat))]

melt_plot_mat <- melt(plot_mat)


p<-ggplot(data=melt_plot_mat,aes(x=Var2, y=Var1, fill=value))+
      geom_tile(color = "gray")+
      scale_fill_gradient(low = "white", high = "blue",limits = c(0,0.2))


print(p)



#Split into high and low mutation samples
select_high <- which(colSums(mutational_matrix) > 2 * median(colSums(mutational_matrix)))
select_low  <- which(colSums(mutational_matrix) < 2 * median(colSums(mutational_matrix)))

mutational_matrix_high <- mutational_matrix[, select_high]
mutational_matrix_low <- mutational_matrix[, select_low]


#Cluster rows based on mutational_matrix
sample_order <- hclust(dist(t(mutational_matrix_high)),method="complete")$order
sample_order <- colnames(mutational_matrix_high)[sample_order]
mutational_matrix_high <- mutational_matrix_high[, match(sample_order ,colnames(mutational_matrix_high))]



#Create a cosine similarity matrix
cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix_high, cosmic_signatures)

#Fit mut_matrix to cosmic
fit_to_cosmic <- fit_to_signatures(mutational_matrix_high,cosmic_signatures)

#Filter signatures to plot
select <- which(rowSums(fit_to_cosmic$contribution)/sum(rowSums(fit_to_cosmic$contribution)) > 0.01)

#Plot
print(plot_contribution(fit_to_cosmic$contribution[select ,],cosmic_signatures[,select],coord_flip= TRUE,mode="absolute"))


#Same for low


#Cluster rows based on mutational_matrix
sample_order <- hclust(dist(t(mutational_matrix_low)),method="complete")$order
sample_order <- colnames(mutational_matrix_low)[sample_order]
mutational_matrix_low <- mutational_matrix_low[, match(sample_order ,colnames(mutational_matrix_low))]



#Create a cosine similarity matrix
cos_sim_samples_cosmic <- cos_sim_matrix(mutational_matrix_low, cosmic_signatures)

#Fit mut_matrix to cosmic
fit_to_cosmic <- fit_to_signatures(mutational_matrix_low,cosmic_signatures)

#Filter signatures to plot
select <- which(rowSums(fit_to_cosmic$contribution)/sum(rowSums(fit_to_cosmic$contribution)) > 0.01)

#Plot

print(plot_contribution(fit_to_cosmic$contribution[select ,],cosmic_signatures[,select],coord_flip= TRUE,mode="absolute"))









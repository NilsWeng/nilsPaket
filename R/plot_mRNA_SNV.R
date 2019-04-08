#' plot_mRNA_SNV
#'
#' Function that given a list of genes and samples plot the CNV/SNV correlation heatmap 
#' @param samples List of sample-id 
#' @param genes List of Hugo- gene id:s
#' @param cluster_names Plot cluster id on left side? False by default
#' @param mRNA_DF Data frame containing mRNA data( from mRNA_exp(genes))
#' @export
#' @examples
#' plot_mRNA_SNV(c(XX,YY),c("ACC","BRCA1"),"Good title",cluster_names=FALSE,mRNA_DF)

#' @import ggplot2
#' @import dplyr
#' @import grid
#' @importFrom plyr join
#' @importFrom reshape2 melt
#' @importFrom zoo rollmean
#' @importFrom data.table setcolorder




plot_mRNA_SNV <- function(samples,genes,cluster_names=FALSE,mRNA_DF){
  
  

  
  setwd(main_wd)
  

  #Load MC3 data
  load("MC3.rda")
  MC3 <- MC3 %>% select("Hugo_Symbol","Chromosome","Start_Position"
                        ,"End_Position","Tumor_Sample_Barcode","Variant_Classification")
  
  #Filter on genes
  MC3 <- MC3[MC3$Hugo_Symbol %in% genes ,]
  
  #Modify sample-id
  MC3$Tumor_Sample_Barcode <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","", MC3$Tumor_Sample_Barcode)
  
  #Filter on samples
  MC3 <- MC3[MC3$Tumor_Sample_Barcode %in% samples ,]
  
  
  #Common genes
  common_genes <- intersect(unique(MC3$Hugo_Symbol),rownames(mRNA_DF))
  common_genes <- common_genes[order(match(common_genes,genes))]
  common_genes_DF <- data.frame("Hugo_Symbol"=common_genes)
  


  
  MC3 <- MC3[MC3$Hugo_Symbol %in% common_genes ,]
  mRNA_DF <- mRNA_DF[rownames(mRNA_DF) %in% common_genes ,]
  



  
  #Select only  samples common between MC3 and mRNA
  common_samples <- intersect(MC3$Tumor_Sample_Barcode,colnames(mRNA_DF))
  #Set order to match samples
  common_samples <- samples[samples %in% common_samples]
  
  
  
  MC3 <- MC3[MC3$Tumor_Sample_Barcode %in% common_samples ,]



  
  #mRNA section
  
  
  #mRNA_DF <- mRNA_DF[order(match(mRNA_DF$`Gene Symbol`,common_genes)) ,]
  mRNA_DF <- as.data.frame(t(mRNA_DF))
  mRNA_DF <- setcolorder(mRNA_DF,common_genes)
  mRNA_DF <- t(mRNA_DF)
  
  
  
  cancer_vector <- get_cancer_type(colnames(mRNA_DF))
  cancerDF <- data.frame("Cancer"=cancer_vector,"sample"=colnames(mRNA_DF))
  
  
  
  
  mRNA_exp.m <-  matrix(nrow=length(common_samples),ncol = length(common_genes))
  colnames(mRNA_exp.m) <- common_genes
  rownames(mRNA_exp.m) <- common_samples

  for (sample in common_samples){
    
    
    cancer_type <- as.character(cancerDF[cancerDF$sample == sample ,1])
    cancer_samples <- cancerDF %>% filter(Cancer==cancer_type) %>% select(sample)
    cancer_samples <- as.character(cancer_samples$sample)
    cancer_samples <- mRNA_DF[, colnames(mRNA_DF) %in% cancer_samples]
    cancer_samples <- as.matrix(cancer_samples)
    mode(cancer_samples) <- "numeric"
    
   
    #Normalization / Standardisation
    norm_cancer_samples <- log2(cancer_samples + 1)
    #standard_cancer_samples <- (norm_cancer_samples - rowMeans(norm_cancer_samples))/(sd(rowMeans(norm_cancer_samples)))
    Z_standard_cancer_samples <- apply(norm_cancer_samples, 1, function(x)(x-mean(x))/(sd(x)))
    
    
    tresholded_cancer <-  Z_standard_cancer_samples
    tresholded_cancer1 <-  Z_standard_cancer_samples
    
    
    
    tresholded_cancer1[which(tresholded_cancer <= -3)] <- "-3S"
    tresholded_cancer1[which(tresholded_cancer >=  3)] <- "+3S"
    tresholded_cancer1[which(tresholded_cancer >=  2 & tresholded_cancer <  3)] <- "+2S"
    tresholded_cancer1[which(tresholded_cancer <= -2 & tresholded_cancer >  -3)] <- "-2S"
    tresholded_cancer1[which(tresholded_cancer >  -2 & tresholded_cancer <  2)] <- "0"
    
    
    tresholded_cancer <- t(tresholded_cancer1)
    tresholded_cancer <- as.vector(tresholded_cancer[, colnames(tresholded_cancer)==sample])
    
    
    
    
    
    mRNA_exp.m[rownames(mRNA_exp.m) == sample] <- tresholded_cancer
    
    
    
    
    
  }
  
  

  mRNA_exp.m <- t(mRNA_exp.m)
  

  
  
  matrix_frame <-  matrix(nrow=length(common_genes),ncol = length(common_samples))
  colnames(matrix_frame) <- common_samples
  rownames(matrix_frame) <- common_genes
  matrix_frame <- as.data.frame(matrix_frame)
  
  snv.m <- matrix_frame
  
  
  
  #Order after gene order (important for getting the right genes)
  MC3 <- MC3[order(match(MC3$Hugo_Symbol,common_genes))]
 
  
  
  
  
  for (sample in common_samples){
    

    
    #Create vector for SNV
    fill_vector_SNV <- MC3 %>% filter(MC3$Tumor_Sample_Barcode == sample) %>% select(Hugo_Symbol,Variant_Classification)
    #To handle genes with several mutations
    fill_vector_SNV <- unique(fill_vector_SNV) #If a gene have several of one type just write that type
    dup <-duplicated(fill_vector_SNV$Hugo_Symbol)
    dup <- fill_vector_SNV[dup ,1]
    fill_vector_SNV[fill_vector_SNV$Hugo_Symbol %in% dup,2] <- "Several"
    fill_vector_SNV <- unique(fill_vector_SNV)
    fill_vector_SNV <- plyr::join(common_genes_DF,fill_vector_SNV,by="Hugo_Symbol")
    
    
    
    
    snv.m[colnames(snv.m) == sample] <- fill_vector_SNV$Variant_Classification
    
    
  }
  
  
 

  
  snv.m <- as.matrix(snv.m)
  
  snv.m <- t(snv.m)
  mRNA_exp.m <- t(mRNA_exp.m)
  mRNA_exp.m[which(mRNA_exp.m == "NaN")] <- "0"
  
  
  
  # If there are to many types bundle together unessesary into other
  non_AA_mod <- c("3'Flank","3'UTR","5'Flank","5'UTR","Intron","Silent","RNA")
  frame_shift <- c("Frame_Shift_Del","Frame_Shift_Ins")
  expression <- c("Splice_Site","Translation_Start_Site")
  inframe    <-c("In_Frame_Del","In_Frame_Ins")
  AA_mod     <- c("Frame_shift","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation")
  
  
  snv.m[snv.m %in% non_AA_mod]  <- "Other"
  snv.m[snv.m %in% frame_shift] <- "Frame_shift"
  snv.m[snv.m %in% expression]  <- "Expression regulation"
  snv.m[snv.m %in% inframe]     <- "In frame del/ins"
  snv.m[snv.m %in% AA_mod]      <- "AA modifying mut"
  
  
  #Prepare data for plotting
  
  
  
  
  
  
  
  #Prepare data for plotting
  
  #if (plot_id ==TRUE){
    
    #test <- ClusterDF %>% filter(gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",TCGA_code) %in% common_samples)
    #test <- as.vector(test$sample)
    #test <- as.character(test)
    #test <- paste(rep("ID: ",length(test)),test,sep="")
    #row.names(snv.m) <- test
    #row.names(mRNA_exp.m) <- test
  #s}            
  
  
  
  
  
  
  #Create dataframe in format for GGplot
  GGdata <- melt(mRNA_exp.m)
  GGdata1 <- melt(snv.m)
  GGdata$mRNA <- as.factor(GGdata$value)
  GGdata$snvs <- as.factor(GGdata1$value)
  GGdata <- GGdata %>% select(Var1,Var2,mRNA,snvs)
  colnames(GGdata) <- c("sample","gene","mRNA","snvs")
  
  
  
  snv_shapes <-  c(0,1,2,3,6)
  textcol <- "grey40"
  
  #In order to handle sample-selection missing one CNV type
  mRNA_types <- c("-3S","-2S","0","+2S","+3S")
  mRNA_col <- c("#E69F00","#D55E00","#56B4E9","#0072B2","grey70")
  mRNA_col_DF <- data.frame("type"=mRNA_types,"colour"=mRNA_col)
  types_in_samples <- unique(names(table(mRNA_exp.m)))
  mRNA_col_DF <- mRNA_col_DF %>% filter(type %in% types_in_samples) %>% select(colour)
  mRNA_col <- as.vector(mRNA_col_DF$colour)
  
  
  #For getting vertical lines marking clusters
  test <- ClusterDF %>% filter(gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",TCGA_code) %in% common_samples)
  
  hline_pos <- c()
  for (cluster in unique(test$cluster)){
    
    pos <- tail(which(test$cluster == cluster),n=1)
    pos <- pos + 0.5
    
    hline_pos <- c(hline_pos,pos)
    
  }
  
  unique(GGdata$mRNA)
  
  
  
  
  
  
  
  P <- ggplot(GGdata,aes(x=gene,y=sample,fill=mRNA))+
    geom_tile()+
    #redrawing tiles to remove cross lines from legend
    geom_tile(colour="white",size=0.25)+
    #remove axis labels, add title
    labs(x="",y="",title="mRNA/SNV")+
    #remove extra space
    #scale_y_discrete(expand=c(0,0),labels = y_names)+  
    scale_y_discrete(expand=c(0,0))+
    #custom breaks on x-axis
    scale_x_discrete(expand=c(0,0))+
    #custom colours for cut levels and na values
    scale_fill_manual(values=mRNA_col)+
    #Add snv data (remove NA(ie no mutation))
    geom_point(data=subset(GGdata,!snvs=="NA"), aes(shape=snvs),size=1)+
    scale_shape_manual(values=snv_shapes)+
    #Add vertical lines where clusters are
    geom_hline(yintercept=hline_pos)+
    theme(
      #remove legend title
      legend.title=element_blank(),
      #remove legend margin
      legend.spacing  = grid::unit(0,"cm"),
      #change legend text properties
      legend.text=element_text(colour=textcol,size=7,face="bold"),
      #set a slim legend
      legend.key.width=grid::unit(0.4,"cm"),
      #set x axis text size and colour
      axis.text.x=element_text(colour="black",angle=90, hjust=1,vjust = 0.5,size=5),
      #set y axis text colour and adjust vertical justification
      axis.text.y=element_text(vjust = 0.5,size=5),
      #change axis ticks thickness
      axis.ticks=element_line(size=0.4),
      #change title font, size, colour and justification
      plot.title=element_text(colour=textcol,hjust=0,size=14,face="bold"),
      plot.margin=unit(c(0.5,0.5,0.5,3),"cm"),
      panel.border = element_rect(fill = NA))
  
  
  print(P)
  
  if (cluster_names==TRUE){
    
    plot.margin=unit(c(0.5,0.5,0.5,3),"cm")
  }
  
  if (cluster_names==TRUE){
    
    library(grid)
    library(zoo)
    top <- 0.96
    bot <- 0.09
    
    units <- (top-bot)/length(common_samples)
    hline_pos1 <- c(0,hline_pos)
    hline_pos1 <- rollmean(hline_pos1,k=2)
    hline_pos1 <- units * hline_pos1
    hline_pos1 <- hline_pos1 + 0.08
    
    
    my_text <- paste(rep("Cluster" , length(unique(ClusterDF$cluster))),unique(ClusterDF$cluster),sep=" ")
    #my_text <- paste(rep("Cluster",14),c(1:14),sep=" ")
    
    grid.text(my_text,x=0.05, y=hline_pos1[1:14],gp=gpar(col="firebrick", fontsize=10, fontface="bold"))
    
    
  }
  
  
  
  
  
}

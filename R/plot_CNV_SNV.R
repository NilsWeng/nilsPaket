#' plot_CNV_SNV
#'
#' Function that given a list of genes and samples plot the CNV/SNV correlation heatmap 
#' @param samples List of sample-id 
#' @param genes List of Hugo- gene id:s
#' @param title1 Plot title
#' @param plot_id FALSE <- Use TCGA code for plotting, TRUE <- use number
#' @param cluster_names Plot cluster id on left side? False by default
#' @export
#' @examples
#' plot_CNV_SNV(c(XX,YY),c("ACC","BRCA1"),"Good title",cluster_names=FALSE)

#' @import ggplot2
#' @import dplyr
#' @import grid
#' @importFrom plyr join
#' @importFrom reshape2 melt
#' @importFrom zoo rollmean




plot_CNV_SNV <- function(samples,genes,title1,plot_id,cluster_names=FALSE){
  
  
  
  
  
  #Read GISTIC-files containing copynumber calls
  
  setwd(main_wd)
  file_list <- list.files(pattern="all_thresholded.by_genes.txt",recursive = TRUE)
  #file_list <- list.files(pattern="all_data_by_genes.txt",recursive = TRUE)
  GISTIC <- read_GISTIC(file_list,genes)
  
  #Treshold GISTIC file
  #GISTIC1 <- GISTIC
  
  #Convert arm significant genes to same as significant
  #GISTIC1[GISTIC1 >= 1 ] <- 1
  #GISTIC1[GISTIC1 <= -1 ] <- -1
  
  #GISTIC1[GISTIC1 >= 1.3 ] <- 2
  #GISTIC1[GISTIC1<1.3 & GISTIC1>-1.1] <- 0
  #GISTIC1[0.7 < GISTIC1 & GISTIC1 < 1.7]<- 1
  #GISTIC1[-0.7 <= GISTIC1 & GISTIC1 <= 0.7]<- 0
  #GISTIC1[-1.3 < GISTIC1 & GISTIC1 < -0.7] <- -1
  #GISTIC1[GISTIC1 <= -1.1] <- -2
  #GISTIC[, 4:ncol(GISTIC)] <-  GISTIC1[, 4:ncol(GISTIC1) ]
  
  
  #select on samples present in samples
  keep <- c(samples,c("Gene Symbol","Locus ID","Cytoband"))
  colnames(GISTIC) <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(GISTIC))
  GISTIC <- GISTIC[, colnames(GISTIC) %in% keep]
  
  
  #Load MC3 data
  load("MC3.rda")
  
  MC3 <- MC3 %>% select("Hugo_Symbol","Chromosome","Start_Position"
                        ,"End_Position","Tumor_Sample_Barcode","Variant_Classification")
  
  #Filter on genes
  MC3 <- MC3[MC3$Hugo_Symbol %in% genes ,]
  #Modify sample-id
  MC3$Tumor_Sample_Barcode <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","", MC3$Tumor_Sample_Barcode)
  
  #Common genes
  common_genes <- intersect(unique(MC3$Hugo_Symbol),GISTIC$`Gene Symbol`)
  
  MC3 <- MC3[MC3$Hugo_Symbol %in% common_genes ,]
  GISTIC <- GISTIC[GISTIC$`Gene Symbol` %in% common_genes ,]
  
  #Filter on samples
  MC3 <- MC3[MC3$Tumor_Sample_Barcode %in% samples ,]
  
  
  #Select only genes and samples common between MC3 and GISITC
  #Common samples
  common_samples <- intersect(MC3$Tumor_Sample_Barcode,colnames(GISTIC))
  
  MC3 <- MC3[MC3$Tumor_Sample_Barcode %in% common_samples ,]
  GISTIC <- GISTIC %>% select("Gene Symbol","Locus ID","Cytoband",common_samples)

  
  #MC3 <- MC3[order(MC3$Hugo_Symbol) ,]
  #GISTIC <- GISTIC[order(GISTIC$`Gene Symbol`) ,]
  
  #Order after order in input vector of samples
  common_samples <- common_samples[order(match(common_samples,samples))]
  
  
  
  
  #Create a matrix of SNVs and CNVs
  common_genes <- common_genes[order(match(common_genes,genes))]
  common_genes_DF <- data.frame("Hugo_Symbol"=common_genes)
  
  
  
  matrix_frame <-  matrix(nrow=length(common_genes),ncol = length(common_samples))
  colnames(matrix_frame) <- common_samples
  rownames(matrix_frame) <- common_genes
  matrix_frame <- as.data.frame(matrix_frame)
  
  cnv.m <- matrix_frame
  snv.m <- matrix_frame
  

  
  #Order after gene order (important for getting the right genes)
  MC3 <- MC3[order(match(MC3$Hugo_Symbol,common_genes))]
  GISTIC <- GISTIC[order(match(GISTIC$`Gene Symbol`,common_genes)) ,]
  
  
  
  
  for (sample in common_samples){
    
    #Create vector for  CN
    fill_vector_CN <- GISTIC %>% select(`Gene Symbol`,sample)
    colnames(fill_vector_CN) <- c("Hugo_Symbol","CN")
    
    #Create vector for SNV
    fill_vector_SNV <- MC3 %>% filter(MC3$Tumor_Sample_Barcode == sample) %>% select(Hugo_Symbol,Variant_Classification)
    #To handle genes with several mutations
    fill_vector_SNV <- unique(fill_vector_SNV) #If a gene have several of one type just write that type
    dup <-duplicated(fill_vector_SNV$Hugo_Symbol)
    dup <- fill_vector_SNV[dup ,1]
    fill_vector_SNV[fill_vector_SNV$Hugo_Symbol %in% dup,2] <- "Several"
    fill_vector_SNV <- unique(fill_vector_SNV)
    fill_vector_SNV <- plyr::join(common_genes_DF,fill_vector_SNV,by="Hugo_Symbol")
    
    
    
    
    cnv.m[colnames(cnv.m) == sample] <- fill_vector_CN$CN
    snv.m[colnames(snv.m) == sample] <- fill_vector_SNV$Variant_Classification
    
    
  }
  
  
  cnv.m <- as.matrix(cnv.m)
  class(cnv.m) <- "numeric"
  cnv.m   <- cnv.m + 2
  
  
  snv.m <- as.matrix(snv.m)
  
  cnv.m <- t(cnv.m)
  snv.m <- t(snv.m)
  
  
  # Catagorize snvs ---------------------------------------------------------
  
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
  
  if (plot_id ==TRUE){
    
    test <- ClusterDF %>% filter(gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",TCGA_code) %in% common_samples)
    test <- as.vector(test$sample)
    test <- as.character(test)
    test <- paste(rep("ID: ",length(test)),test,sep="")
    row.names(snv.m) <- test
    row.names(cnv.m) <- test
  }            
    
 
  
  
  
  
  #Create dataframe in format for GGplot
  GGdata <- melt(cnv.m)
  GGdata1 <- melt(snv.m)
  GGdata$cnvs <- as.factor(GGdata$value)
  GGdata$snvs <- as.factor(GGdata1$value)
  GGdata$cnvs <- as.factor(GGdata$value)
  GGdata <- GGdata %>% select(Var1,Var2,cnvs,snvs)
  


  colnames(GGdata) <- c("sample","gene","cnvs","snvs")
  

  
  snv_shapes <-  c(0,1,2,3,6)
  textcol <- "grey40"
  
  #In order to handle sample-selection missing one CNV type
  cnv_types <- c(0,1,2,3,4)
  cnv_col <- c("#D55E00", "#E69F00","grey70","#56B4E9","#0072B2")
  cnv_col_DF <- data.frame("type"=cnv_types,"colour"=cnv_col)
  types_in_samples <- sort(unique(c(cnv.m)),decreasing=FALSE)
  cnv_col_DF <- cnv_col_DF %>% filter(type %in% types_in_samples) %>% select(colour)
  cnv_col <- as.vector(cnv_col_DF$colour)
  
  
  #For getting vertical lines marking clusters
  test <- ClusterDF %>% filter(gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",TCGA_code) %in% common_samples)
  
  hline_pos <- c()
  for (cluster in unique(test$cluster)){
    
    pos <- tail(which(test$cluster == cluster),n=1)
    pos <- pos + 0.5
    
    hline_pos <- c(hline_pos,pos)
    
  }
  
  #Cancer id as y axis labels
  if(FALSE){
    samples_id <- data.frame("sample"=rownames(cnv.m))
    test <- ClusterDF %>% select(TCGA_code,Study.Abbreviation)
    colnames(test) <- c("sample","cancer")
    test$sample <- gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",test$sample)
    samples_id <- left_join(samples_id,test,by="sample")
    
  }

  
  
  
  P <- ggplot(GGdata,aes(x=gene,y=sample,fill=cnvs))+
    geom_tile()+
    #redrawing tiles to remove cross lines from legend
    geom_tile(colour="white",size=0.25)+
    #remove axis labels, add title
    labs(x="",y="",title=title1)+
    #remove extra space
    scale_y_discrete(expand=c(0,0))+  
    #scale_y_discrete(labels=samples_id$cancer)+
    #custom breaks on x-axis
    scale_x_discrete(expand=c(0,0))+
    #custom colours for cut levels and na values
    scale_fill_manual(values=cnv_col)+
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

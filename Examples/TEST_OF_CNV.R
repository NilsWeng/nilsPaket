#Example pipeline
rm(list=ls())
gc()



#Set wd to folder containing all data
main_wd <- "C:/Users/Nils_/OneDrive/Skrivbord/Main/Data"
setwd(main_wd)

file_list <- list.files(pattern="all_thresholded.by_genes.txt",recursive = TRUE)

genes <- c("FGFR1","PIK3CA")
GISTIC <- read_GISTIC(file_list,genes)

cancer_vector <- get_cancer_type(gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",colnames(GISTIC[,4:ncol(GISTIC)])))

luad <- colnames(GISTIC)[which(cancer_vector %in% c("LUAD"))]
lusc <- colnames(GISTIC)[which(cancer_vector %in% c("LUSC"))]




luad_GISTIC <- GISTIC %>% select(luad)
luad_GISTIC <- t(luad_GISTIC)
colnames(luad_GISTIC) <- GISTIC$`Gene Symbol`
lusc_GISTIC <- GISTIC %>% select(lusc)
lusc_GISTIC <- t(lusc_GISTIC)
colnames(lusc_GISTIC) <-GISTIC$`Gene Symbol`
lusc_GISTIC <- as.data.frame(lusc_GISTIC)
luad_GISTIC <- as.data.frame(luad_GISTIC)




count_luad_FGF <- luad_GISTIC %>% dplyr::count(FGFR1)
count_luad_FGF$gene <- rep("FGFR1",nrow(count_luad_FGF))


count_luad_PIK <- luad_GISTIC %>% dplyr::count(PIK3CA)
count_luad_PIK$gene <- rep("PIK",nrow(count_luad_PIK))

colnames(count_luad_PIK ) <- c("type","count","gene")
colnames(count_luad_FGF ) <- c("type","count","gene")


LUAD_DF <- rbind(count_luad_FGF,count_luad_PIK)

LUAD_DF <- LUAD_DF %>% arrange(type)

library(ggplot2)

ggplot(data=LUAD_DF, aes(x=gene, y=count, fill=as.character(type))) +
  geom_bar(stat="identity")





count_lusc_FGF <- lusc_GISTIC %>% dplyr::count(FGFR1)
count_lusc_FGF$gene <- rep("FGFR1",nrow(count_lusc_FGF))


count_lusc_PIK <- lusc_GISTIC %>% dplyr::count(PIK3CA)
count_lusc_PIK$gene <- rep("PIK",nrow(count_lusc_PIK))

colnames(count_lusc_PIK ) <- c("type","count","gene")
colnames(count_lusc_FGF ) <- c("type","count","gene")


lusc_DF <- rbind(count_lusc_FGF,count_lusc_PIK)



library(ggplot2)

ggplot(data=lusc_DF, aes(x=gene, y=count, fill=as.character(type))) +
  geom_bar(stat="identity")





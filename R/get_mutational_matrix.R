#' get_mutational_matrix
#'
#' Function that given vector of samples and genes and a mRNA expression table(from mRNA_exp function) returns a matrix with exp/avg exp
#' @param VCF_files List of directories to VCF files to extract a mutational matrix from
#' @param file_name Name of file you want to createBothm mutmat and GRange object.
#' @export
#' @examples
#' get_mutational_matrix(VCF_files,"Run3")
#' @import MutationalPatterns
#' @import dplyr 






get_mutational_matrix <- function(VCF_files,file_name){
  
  setwd(paste(main_wd,"/VCF_files",sep=""))
  
  #Convert VCF to GRange object
  print("Converting VCF to GRange-object")
  
  name_vector  <- gsub("/.*","",VCF_files)
  GRange_vcf <- read_vcfs_as_granges(VCF_files,name_vector,ref_genome)
  
  
  
  #Extract mutational matrix
  print("Creating mutational matrix")
  mutational_matrix <-  mut_matrix(GRange_vcf,ref_genome)
  
  setwd(paste(main_wd,"/mut_matrix",sep=""))
  save(mutational_matrix,file=paste("mut_matrix_",file_name,".rda",sep=""))
  save(GRange_vcf,file=paste("GRange_",file_name,".rda",sep=""))
}






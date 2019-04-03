#' get_cancer_type
#'
#' Function that given a vector of samples return a vector of their cancer type
#' @param samples Samples one wants to find cancer type of.
#' @export
#' @examples
#' get_cancer_type("TCGA-BT-A2LB-01A" "TCGA-DK-A1AC-01A" "TCGA-DK-A3WW-01A")

#' @importFrom dplyr left_join
#' @return TSS_code





get_cancer_type <- function(TCGA_code){
  
  current_wd <- getwd()
  setwd(main_wd)
  TSS2Study <- read.table("TSS2Studyabb.txt",header=TRUE,stringsAsFactors = FALSE)
  #Handle NA -> USCS , read as na value
  TSS2Study[which(is.na(TSS2Study$TSS.Code)),1] <- "NA"
  
  
  TSS_code <- gsub("TCGA-","",gsub("-[A-Z0-9]*-[A-Z0-9]*$","",TCGA_code))
  TSS_code <- left_join(data.frame("TSS.Code"=TSS_code),TSS2Study,by="TSS.Code")
  TSS_code <- as.vector(TSS_code$Study.Abbreviation)
  
  setwd(current_wd)
  return(TSS_code)
  
  
}
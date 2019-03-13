#' MC3_to_VCF
#'
#' MC3_to_VCF
#'
#'
#' @export
#' @examples
#' MC3_to_VCF()

#' @import dplyr
#' @import data.table
#' @import plyr



MC3_to_VCF <- function(){
  
  #Lookup table for TSS code to cancer type
  TSS2Study_DF <- read.table("TSS2Studyabb.txt",header=TRUE)
  
  MC3_DF <- fread('mc3.v0.2.8.PUBLIC.maf.gz')
  MC3_DF <- select(MC3_DF, Chromosome,Tumor_Seq_Allele2,Tumor_Sample_Barcode,Start_Position,Reference_Allele)
  MC3_DF$TSS.Code <-gsub("-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*-[A-Z0-9]*$","",gsub("TCGA-","",MC3_DF$Tumor_Sample_Barcode))
  MC3_DF <- join(MC3_DF,TSS2Study_DF, by="TSS.Code")
  
  #Save mc3-file in separate folders as VCF-files
  
  MC3_DF_CancerType <- split (MC3_DF, f=MC3_DF$Study.Abbreviation,drop=TRUE)
  
  
  vcffile=vector()
  subtype=vector()
  sample=vector()
  
  
  
  #Start of loops --------------------------------------------
  
  for (cancerType in MC3_DF_CancerType){ # Loop over all cancer types
    
    #Add a catch that prevents closing loop
    if (dim(cancerType)[1]  == 0) {
      
      next()
      
    }
    
    
    
    directory_name <- paste(dir_name,"/",unique(cancerType$Study.Abbreviation),sep="") 
    
    #crate folder if it doesnt exist
    ifelse(!dir.exists(file.path(directory_name)), dir.create(file.path(directory_name)), FALSE)
    print (directory_name)
    setwd(directory_name)
    
    
    #Split cancer type into individual samples
    MC3_DF_Sample <- split (cancerType, f= cancerType$Tumor_Sample_Barcode,drop = TRUE) # not sure if this works
    
    
    
    
    
    
    
    for (Sample in MC3_DF_Sample){ #Loop over each sample in respective cancer type
      
      
      if (dim(Sample)[1]  == 0) { ########### Why does CancerType get split into empty samples??????????!!!!!!!!
        #Seem to split MC3_DF_Sample but still save old levels  # try drop = TRUE in split command
        
        next()
        
      } 
      
      
      
      sample_id <- unique(Sample$Tumor_Sample_Barcode);
      vcfdata <- matrix(".", nrow=nrow(Sample), ncol = 10);
      columns=c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",sample_id )
      colnames(vcfdata) <- columns
      
      vcfdata[,1]=as.character(Sample$Chromosome)
      vcfdata[,2]=as.character(Sample$Start_Position)
      vcfdata[,4]=as.character(Sample$Reference_Allele)
      vcfdata[,5]=as.character(Sample$Tumor_Seq_Allele2)
      vcfdata[,9]="GT"
      vcfdata[,10]="1/0"
      #vcfdata[,3]=as.character(sample_id)
      outfile <- paste(sample_id,".vcf", sep="")
      #outfile = paste(directory_name,outfile, sep="/")
      write.table(vcfdata, file=outfile, row.names=FALSE, sep="\t", quote=FALSE)
      
      subtype <- c(subtype,as.character(unique(cancerType$Study.Abbreviation)))
      sample <- c(sample,as.character(sample_id))
      vcffile <- c(vcffile, paste (unique(cancerType$Study.Abbreviation),outfile,sep = "/") )
      
    }
    
    
  }
  
  
  
  cohort <- data.frame("vcf"=vcffile, "subtype"=subtype, "sample"=sample) 
  setwd(main_wd)
  write.table(cohort,file="Cohort.txt")
  
  
  
}
  
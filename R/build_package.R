#' Build package
#'
#' Build package
#'
#' @param package_name Package to document,install and load
#'
#' @return None
#' @export
#' @examples
#' build_package(nilsPaket)

#' @import devtools
#' @import roxygen2


build_package <- function(){
    
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/GITrepo/nilsPaket")
  document()
  setwd("C:/Users/Nils_/OneDrive/Skrivbord/GITrepo")
  install("nilsPaket")
  library(nilsPaket)
  print(done)
  
 
  
}


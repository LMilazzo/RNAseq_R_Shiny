
AssertFile <- function(file){
  
  if(is.null(file)){
    return(FALSE)
  }
  
  if(grepl(".tsv", file)){
    r <- read.csv(file, sep = "\t")
    return(r)
  }
  
  if(grepl(".csv", file)){
    r <-  read.csv(file)
    return(r)
  }

  return(FALSE)
  
}
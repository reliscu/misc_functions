upper_first <- function(string){
  
  split <- strsplit(string, "")[[1]]
  return(paste0(toupper(split[1]), paste(split[2:length(split)], collapse="")))
  
}
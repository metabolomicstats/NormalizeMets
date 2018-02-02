reducedlabel <- function(namevec, maxnames=40) {
  #check if namevec is longer than maxnames
  
  if (length(namevec) <= maxnames){
    return(namevec)
  }
  
  #create blank vector to host new names
  outnames <- c(rep("",length(namevec)))
  
  #find how often to show real name
  skiplen <- ceiling(length(namevec)/maxnames)
  
  for (i in 1:length(namevec)){
    
    if((i%%skiplen)==0){
      outnames[i] <- namevec[i]
    }
  }
  
  return(outnames)
  
  
}
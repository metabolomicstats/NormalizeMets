#makes a scale of round numbers given min, max and the number of steps
pretty_simp <- function(min,max,n=8){

  max_n <- max
  min_n <- min
  factor <- 1
  
  while((max_n - min) >= 10^6){
    factor <- factor*10
    max_n <- ceiling(max_n/10)
    min_n <- floor(min_n/10)
  }
  
  return(factor*(pretty(min_n:max_n,n)))
  
}
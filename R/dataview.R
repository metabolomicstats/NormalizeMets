#' Data Viewing 
#' 
#' Prints out a trimmed version of the data matrix
#' 
#' @param data a data frame or a data table
#' 
#' @return The trimmed data
#' 
#' @author Gavriel Olshanksy and Alysha M De Livera
#' 
#' @examples 
#' data("alldata_eg")
#' featuredata_eg<-alldata_eg$featuredata
#' dataview(featuredata_eg)
#' 
# A function to show some of the data for basic understanding/viewing of the data
#
#'@export dataview
dataview <- function(data){
  if (class(data) %in% c("list")) 
    stop("Input data must be a vector, a matrix, or a dataframe")
  data <- as.data.frame(data)
  hdata <- data
  dsize <- dim(hdata)
  
  # Check and update rows 
  if (dsize[1] > 10){
    hdata <- hdata[1:10,]
    
  }
  # Check and update columns
  if (dsize[2] > 8){
    hdata <- hdata[,1:9]
    hdata[,9] <- rep("...",nrow(hdata))
    colnames(hdata)[9] <- "..."
  }
  write(
    paste('\n - The data consists of ',dsize[1],' rows by ',dsize[2],' columns \n',sep = "")
    ,
    ''
  )
  return(hdata)
}
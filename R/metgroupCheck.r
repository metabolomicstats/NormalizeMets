#' Check data compatibility
#' 
#' Check if the group data matches the metabolite data in size, returns error
#' message if not compatible
#' 
#' 
#' @param metdata A data frame in the met data format. This should have sample
#' names in the first column to be read as row names and the variables in the
#' remaining columns. These variables can be metabolites, masses, retention
#' times, bins, areas or any other metabolomics variables of interest.
#' @param groupdata A data frame or a table with optional sample names in the 
#' first column to be read as row names and group names in the following column.
#' @author Alysha M De Livera, Gavriel Olshansky
#' 
#'
#' 
#' 
metgroupCheck <- function(metdata, groupdata)
{
 #  make sure equal number of rows in both metdata and gorupdata
  rowsmet <- nrow(metdata)
  rowsgroup <- length(groupdata)
  if (rowsmet != rowsgroup){
    stop(
      paste("metdata needs to have same number of rows as groupdata")
    )
  }
  return()
}

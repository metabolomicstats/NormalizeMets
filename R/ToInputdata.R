#' Convert to Inputdata
#' 
#' Returns an inputdata type data frame by combining group information from
#' groupdata with metabolites information from featuredata
#' 
#' 
#' @param featuredata featuredata A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' See NormalizeMets Vignette for details.
#' @param groupdata A data frame or a vector with group names.
#' @return \item{inputdata}{A data frame in the input data format. This has
#' sample names in the first column, group names in the second column, and the
#' metabolite variables in the remaining columns. }
#' @author Alysha M De Livera, Gavriel Olshansky
ToInputdata <- function(featuredata,groupdata)
{
  
  # get groups from groupdata
  groups <- as.data.frame(groupdata)
  
  # Check if Data compatible
  metgroupCheck(featuredata,groupdata)
  
  # make inputdata format data frame
  inputdata <- cbind(groups,featuredata)
  
  return(inputdata)
}

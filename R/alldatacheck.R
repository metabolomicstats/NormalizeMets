#' Check all data
#' 
#'Some preliminary checks of the data formats 
#'for the NormalizeMets package 
#' 
#' @param featuredata featuredata A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' See NormalizeMets Vignette for details.
#' @param sampledata A dataframe that contains sample specific information. 
#' Unique sample names should be provided as row names. 
#' @param metabolitedata A dataframe that contains metabolite specific information. 
#' @author Alysha M De Livera, Gavriel Olshansky
#alldatacheck function
alldatacheck<-function(featuredata,
                       sampledata=NULL, 
                       metabolitedata=NULL){
  
  #local function to check if names are unique
  is.unique <- function(my.names){
    return(all(!duplicated(my.names)))
  }
  
  #Check featuredata    
  if (!class(featuredata) %in% c("data.frame", "matrix")) 
    stop("featuredata must be a matrix or a dataframe")
  if (is.null(rownames(featuredata)))
    stop("featuredata must have row names") 
  else if (!is.unique(rownames(featuredata)))
    stop("featuredata rownames must be unique")
  if (is.null(colnames(featuredata)))
    stop("featuredata must have column names")
  else if (!is.unique(colnames(featuredata)))
    stop("featuredata column names must be unique")
  
  #Check sampledata
  if (!is.null(sampledata)){
    if (class(sampledata) != c("data.frame")) 
      stop("sampledata must be a dataframe")
    if (is.null(rownames(sampledata)))
      stop("sampledata must have row names") 
    if (!identical(rownames(sampledata),rownames(featuredata)))
      stop("featuredata and sampledata row names must be identical")
  }
  
  if (!is.null(metabolitedata)){
    if (class(metabolitedata) != c("data.frame")) 
      stop("metabolitedata must be a dataframe")
    if (is.null(rownames(metabolitedata)))
      stop("sampledata must have row names") 
    if (!identical(rownames(metabolitedata),colnames(featuredata)))
      stop("featuredata column names and metabolitedata row names must be identical")
  }
}

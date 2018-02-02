#' Log transformation
#' 
#' Log transform a metabolomics feature data matrix.
#' 
#' 
#' @param featuredata featuredata A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' See NormalizeMets Vignette for details.
#' @param base The base with respect to which logarithms are computed. The
#' default computes the natural logarithm.
#' @param zerotona A logical indicating whether any zeros should be converted to missing prior to
#' log transforming. By default, this is set to FALSE.
#' @param saveoutput A logical indicating whether the output should be saved as
#' a .csv file.
#' @param outputname The name of the output file if the output has to be saved.
#' @return The result is an object of class
#' \code{\link[NormalizeMets:alldata]{alldata}}. 
#' @author Alysha M De Livera
#' @examples
#' 
#'     data(alldata_eg)
#'     lg <- LogTransform(alldata_eg$featuredata)
#'     dataview(lg$featuredata)
#' 
#' @export LogTransform
LogTransform <- function(featuredata, base=exp(1), saveoutput=FALSE,
    outputname="log.results",zerotona=FALSE)
{
    
    #    Log transform the data for output
    inputdata<-as.matrix(featuredata)
    if (zerotona==TRUE)
      inputdata[which(inputdata==0)]<-NA


    outdata <- log(inputdata, base)
  
    if (saveoutput) {
        write.csv(outdata, paste(c(outputname, ".csv"), collapse=""))
    }
    
    output <- list()
    output$featuredata <- outdata
    output$sampledata <- NULL
    output$metabolitedata <- NULL
    return(structure(output, class="alldata"))
}

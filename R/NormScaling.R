#' Normalisation methods based on scaling
#' 
#' Normalise a metabolomic data matrix according to a specified scaling method.
#' 
#' The normalisation methods based on scaling include normalisation to a total
#' sum, or by the median or mean of each sample, and are denoted by
#' "\code{sum}", "\code{median}", and "\code{mean}" respectively. The method
#' "\code{ref}" normalises the metabolite abundances to a specified reference
#' vector.
#' 
#' @param featuredata featuredata A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' @param method A character string indicating the required scaling-based normalization
#' method. Must be one of "\code{median}", "\code{mean}", "\code{sum}", or
#' "\code{ref}". See NormalizeMets Vignette for details.
#' @param refvec A reference vector to be used with the method "\code{ref}".
#' @param saveoutput A logical indicating whether the normalised data matrix
#' should be saved as a .csv file.
#' @param outputname The name of the output file if the output has to be saved.
#' @param ... Arguments to other functions
#' @return The result is an object of class
#' \code{\link[NormalizeMets:alldata]{alldata}}.
#' @author Alysha M De Livera, Gavriel Olshansky
#' @seealso \code{\link[crmn]{normFit}}.
#' @references De Livera, Alysha M De, M. Aho-Sysi, Laurent Jacob, 
#' J. Gagnon-Bartch, Sandra Castillo, J.A. Simpson, and Terence P. Speed. 
#' 2015. Statistical Methods for Handling Unwanted Variation in 
#' Metabolomics Data. Analytical Chemistry 87 (7). American Chemical Society: 
#' 3606-3615. 
#' 
#' @examples
#' 
#'    ## Reading the data
#'      data(mixdata)
#'      featuredata <- mixdata$featuredata
#'      sampledata<-mixdata$sampledata
#'      metabolitedata<-mixdata$metabolitedata
#'      refvec<-featuredata[,which(metabolitedata$type =="IS")[1]]
#'     
#'     ## Normalise by the median
#'     norm_med <- NormScaling(featuredata, method = "median")
#'     
#'     ## Normalise by a reference vector, in this case an internal standard
#'     norm_is <- NormScaling(featuredata, method = "ref", 
#'         refvec=refvec)
#'
#'     ## Normalise by the sum
#'     norm_sum <- NormScaling(featuredata, method = "sum")
#'
#'    ## Rla Plots after normalising by the median
#'     RlaPlots(norm_med$featuredata, group= sampledata$batch)
#' 
#' @export NormScaling
NormScaling <-function(featuredata, 
                     method=c("median", "mean", "sum", "ref"),
                     refvec=NULL, saveoutput=FALSE, 
                     outputname=NULL,...)
{
  # Match the method
  method <- match.arg(method)
  
  # If method is not one of the above listed, then stop
  if (!is.element(method, 
                  c("median", "mean", "sum", "ref"))
  ) {
    stop("Invalid scaling normalization method")
  }
  
  # Reference vector should be a vector
  if (!is.null(refvec)) {
    if (class(refvec) %in% c("data.frame", "list", "matrix")) {
      stop("Reference should be a vector")
    }
  }
  
  # If there is no refvec is given, get them to enter the reference vector
  if ((method == "ref") & is.null(refvec)) {
    stop("Please enter the reference vector")
  }
  

  # Rename data for processing
  pre_norm <- featuredata
  
  # Median vector
  if (method == "median") {
    norm_vector <- apply(pre_norm, 1, median, na.rm=TRUE)
    # Mean vector
  } else if (method == "mean") {
    norm_vector <- rowMeans(pre_norm, na.rm=TRUE)
    # Sum vector
  } else if (method == "sum") {
    norm_vector <- rowSums(pre_norm, na.rm=TRUE)
  } else if (method == "ref") {
    norm_vector <- refvec
  }
  
  # Subtract the norm_vector
  norm_data <- sweep(pre_norm, 1, norm_vector, "-")
  rownames(norm_data) <- rownames(pre_norm)
  colnames(norm_data) <- colnames(pre_norm)
  
  # make dataframe for output
  outdata <- data.frame(norm_data)
  #Edit column names
  outdata <- editcolnames(outdata)

    
  # Generate the output matrix in .csv format
  if (saveoutput) {
    write.csv(outdata,
              if (!is.null(outputname)) {
                paste(c(outputname, ".csv"), collapse="")
              } else {
                paste(c("normalized_", method, ".csv"), collapse="")
              }
    )
  }


    output <- list()
    output$featuredata <- outdata
    output$sampledata <- NULL
    output$metabolitedata <- NULL
    return(structure(output, class="alldata"))
  }

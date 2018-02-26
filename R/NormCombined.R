#' Normalisation methods based on a combination of methods 
#' 
#' Normalise a metabolomic data matrix using a combination of methods
#' 
#' @param featuredata featuredata A data frame in the featuredata format. 
#' This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' @param methods A character vector indicating which two methods should be used in order. 
#' @param ... Inputs for the \code{NormScaling}, \code{NormQcmets}, and \code{NormQcsamples} 
#' as appropriate. 
#' @param savefinaloutput A logical indicating whether the final normalised data matrix
#' should be saved as a .csv file.
#' @param finaloutputname The name of the final output file if the output has to be saved.
#' @return The result is an object of class
#' \code{\link[NormalizeMets:alldata]{alldata}}. 
#'  
#' @author Alysha M De Livera, Gavriel Olshansky
#' @references De Livera, Alysha M De, M. Aho-Sysi, Laurent Jacob, 
#' J. Gagnon-Bartch, Sandra Castillo, J.A. Simpson, and Terence P. Speed. 
#' 2015. Statistical Methods for Handling Unwanted Variation in 
#' Metabolomics Data. \emph{Analytical Chemistry} 87 (7). American Chemical Society: 
#' 3606-3615.
#' 
#' De Livera, A. M., Dias, D. A., De Souza, D., Rupasinghe, T.,
#' Pyke, J., Tull, D., Roessner, U., McConville, M., Speed, T. P. (2012a)
#' Normalising and integrating metabolomics data. \emph{Analytical Chemistry}
#' 84(24): 1076-10776.
#' 
#' De Livera, A.M., Olshansky, M., Speed, T. P. (2013) Statistical analysis of
#' metabolomics data. \emph{Methods in Molecular Biology} In press.
#' 
#' Gagnon-Bartsch, Johann A., Speed, T. P. (2012) Using control genes to
#' correct for unwanted variation in microarray data. \emph{Biostatistics}
#' 13(3): 539-552.
#' 
#' Redestig, H., Fukushima, A., Stenlund, H., Moritz, T., Arita, M., Saito, K.,
#' Kusano, M. (2009) Compensation for systematic cross-contribution improves
#' normalization of mass spectrometry based metabolomics data. \emph{Analytical
#' Chemistry} 81(19): 7974-7980.
#' 
#' Sysi-Aho, M., Katajamaa, M., Yetukuri, L., Oresic, M. (2007) Normalization
#' method for metabolomics data using optimal selection of multiple internal
#' standards. \emph{BMC Bioinformatics} 8(1): 93.
#' 
#' Dunn,W.B., Broadhurst,D., Begley,P., Zelena,E., Francis-McIntyre,S., Anderson,N., Brown,M., 
#' Knowles,J.D., Halsall,A., Haselden,J.N. et al. (2011) Procedures for 
#' large-scale metabolic profiling of serum and plasma using gas chromatography and 
#' liquid chromatography coupled to mass spectrometry. Nat. Protoc., 6, 1060-1083
#' 
#' @examples
#'  ##Reading the data
#'  data(Didata)
#'  NormCombined(featuredata=Didata$featuredata[order(Didata$sampledata$order),],
#'             sampledata=Didata$sampledata[order(Didata$sampledata$order),],
#'             methods=c("rlsc", "median"),
#'             savefinaloutput=FALSE, 
#'             finaloutputname=NULL)
#'
#' 
#'     
#' @export NormCombined
NormCombined <-function(featuredata,  
                         methods=c("rlsc", "median"),
                         savefinaloutput=FALSE, 
                         finaloutputname=NULL, ...)
{
  
  # If methods are not in the listed, then stop
  if (!all(is.element(methods, 
                      c("is", "nomis", "ccmn", "ruv2",
                        "ruvrand","ruvrandclust",
                        "median", "mean", "sum", "ref", "rlsc")))
  ) {
    stop("Invalid normalization method")
  }
  
  
  # Match the methods
  method1 <- methods[1]
  method2 <- methods[2]
  
  if (method1 %in% c("is", "nomis", "ccmn", "ruv2",
                     "ruvrand","ruvrandclust")){
    out1<-NormQcmets(featuredata, method = method1, ...)
  } else if (method1 %in% c("rlsc")){
    # if (is.null(sampledata)) THESE WARNINGS IN INDIVIDUAL METHODS?
    #   stop("Enter sampledata for the qcrlsc method.")
    out1<-NormQcsamples(featuredata, method = method1,lg=TRUE,...) # sampledata=sampledata, ...)
  } else if (method1 %in% c("median", "mean", "sum", "ref")){
    out1 <-NormScaling(featuredata, method = method1, ...)
  }
  
  
  if (method2 %in% c("is", "nomis", "ccmn", "ruv2",
                     "ruvrand","ruvrandclust")){
    out2<-NormQcmets(out1$featuredata, method = method2, ...)
  } else if (method2 %in% c("rlsc")){
    # if (is.null(sampledata)) THESE WARNINGS IN INDIVIDUAL METHODS?
    stop("rlsc should be done first.")
    #    out2<-NormQcsamples(out1$featuredata, method = method2,...) # sampledata=sampledata, ...)
  } else if (method2 %in% c("median", "mean", "sum", "ref")){
    out2 <-NormScaling(out1$featuredata, method = method2, ...)
  }
  
  # Generate the output matrix in .csv format
  if (savefinaloutput) {
    write.csv(out2$featuredata,
              if (!is.null(finaloutputname)) {
                paste(c(finaloutputname, ".csv"), collapse="")
              } else {
                paste(c("normalized_", method1,"_",method2, ".csv"), collapse="")
              }
    )
  }
  
  if(is.null(out2$sampledata) & !is.null(out1$sampledata))
    out2$sampledata<-out1$sampledata
  
  return(structure(out2, class="alldata"))
  
}


#' Normalisation methods based on quality control metabolites
#' 
#' Normalise a metabolomic data matrix using internal, external standards 
#' and other quality control metabolites
#' 
#' These normalisation methods include "\code{is}" which uses a 
#' single standard, Cross-contribution Compensating
#' Multiple internal standard Normalisation, "\code{ccmn}" (Redestig \emph{et
#' al}., 2009); normalization using optimal selection of multiple internal
#' standards, "\code{nomis}" (Sysi-aho \emph{et al}. 2007), "\code{ruv2}"
#' (De Livera \emph{et al}.  2012a), and "\code{ruvrand}", "\code{ruvrandclust}" 
#' (De Livera \emph{et al}.  2015). 
#' 
#' An overview of these normalisation methods are given by De Livera \emph{et
#' al}. (2015). 
#' 
#' @param featuredata featuredata A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' @param factors For the ccmn method. A vector describing biological factors. 
#' @param factormat For the ruv2 method. A design matrix for the linear model, consisting of biological factors of interest.
#' @param method A character string indicating the required normalization
#' method. Must be one of "\code{is}", "\code{nomis}", "\code{ccmn}", "\code{ruv2}", 
#' "\code{ruvrand} " or "\code{ruvrandclust}". See Details for information.
#' @param isvec A vector of internal standards to be used with the method
#' "\code{is}".
#' @param ncomp Number of PCA components to be used for the "\code{ccmn}"
#' method. If \code{NULL}, this will be determined by cross validation as
#' described by Redestig (2012).
#' @param k Number of factors of unwanted variation to be included in the
#' "\code{ruv}" models.
#' @param plotk For the "\code{ruvrand}" method. A logical indicating whether a bargraph 
#' for the proportion of variance explained by the factors of unwanted variation needs to be plotted
#' @param lambda The regularization parameter for the "\code{ruvrand}" method which depends on k. 
#' If not entered, it will be estimated. See DeLivera et al, 2015 for details.
#' @param qcmets A vector indicating which metabolites should be used as the
#' internal, external standards or other quality control metabolites 
#' in the "\code{ruv}" models, or as multiple internal
#' standards in the "\code{ccmn}" and "\code{nomis}" methods.
#' @param maxIter For the "\code{ruvrandclust}" method. Maximum number of 
#' iterations for "\code{ruvrandclust}" method.
#' @param nUpdate For the "\code{ruvrandclust}" method. Update the unwanted 
#' variation component every nUpdate iterations.
#' @param lambdaUpdate For the "\code{ruvrandclust}" method. A logical indicating whether 
#' the regularization parameter needs to be updated
#' @param p For the "\code{ruvrandclust}" method. The number of clusters to be used 
#' in the k-means clustering.
#'@param saveoutput A logical indicating whether the normalised data matrix
#' should be saved as a .csv file.
#' @param outputname The name of the output file if the output has to be saved.
#' @param ... Other arguments to be passed onto \code{\link[NormalizeMets:LinearModelFit]{LinearModelFit}}.
#' @return If the method is `ruv2', the function will return an object of class 
#' \code{\link[limma:marraylm]{MArrayLM}}, containing F statistics, t statistics, 
#' corresponding confidence intervals, and adjusted and unadjusted p-values. See 
#' \code{\link[NormalizeMets:LinearModelFit]{LinearModelFit}}. For all other methods, the result is an object of class 
#' \code{\link[NormalizeMets:alldata]{alldata}}. Additionally, the list also 
#' contains the removed unwanted variation component (UVcomp),and the
#' results from the optimization algorithm (opt) for the "\code{ruvrandclust}" method
#'  @seealso \code{\link[crmn]{normFit}}.
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
#' standards. \emph{BMC Bioinformatics} 8(1); 93.
#' @examples
#' 
#'    ## Reading the data
#'      data(mixdata)
#'      featuredata <- mixdata$featuredata
#'      sampledata<-mixdata$sampledata
#'      metabolitedata<-mixdata$metabolitedata
#'      isvec<-featuredata[,which(metabolitedata$type =="IS")[1]]
#'      factors<-sampledata$type
#'      qcmets<-which(metabolitedata$type =="IS")
#'     
#'     ## Normalise by an internal or an internal standard
#'     norm_is <- NormQcmets(featuredata, method = "is", isvec=isvec)
#'     PcaPlots(norm_is$featuredata, factors)
#'     
#'     ## Normalise by the NOMIS method
#'     norm_nomis <- NormQcmets(featuredata, method = "nomis", qcmets = qcmets)
#'     PcaPlots(norm_nomis$featuredata, factors)
#'     
#'     ## Normalise by the CCMN method
#'     norm_ccmn <- NormQcmets(featuredata, factors, method = "ccmn", qcmets = qcmets, ncomp = 2)
#'     PcaPlots(norm_ccmn$featuredata, factors)
#'     
#'     ## Normalise using RUV-random method
#'     norm_ruvrand <- NormQcmets(featuredata, method = "ruvrand", qcmets = qcmets, k = 2)
#'     PcaPlots(norm_ruvrand$featuredata, factors)
#'     PcaPlots(norm_ruvrand$uvdata, sampledata$batch, main = "Unwanted batch variation")
#'    
#'     ## Normalise using RUV-random clustering method
#'     ##Not run
#'     #norm_ruvrandclust <- NormQcmets(featuredata, method = "ruvrandclust", qcmets = qcmets, k = 2)
#'     #PcaPlots(norm_ruvrandclust$featuredata, factors)
#'     #PcaPlots(norm_ruvrandclust$uvdata, sampledata$batch, main = "Unwanted batch variation")
#'         
#'     
#' @export NormQcmets
NormQcmets <-function(featuredata, factors = NULL, factormat=NULL,
                      method=c("is", "nomis", "ccmn", "ruv2",
                               "ruvrand","ruvrandclust"),
                      isvec=NULL, ncomp=NULL, k=NULL, plotk=TRUE, 
                      lambda=NULL, qcmets=NULL,  
                      maxIter=200, nUpdate=100,lambdaUpdate=TRUE,
                      p=2,saveoutput=FALSE, 
                      outputname=NULL,...)
{
  # Match the method
  method <- match.arg(method)
  
  # If method is not one of the above listed, then stop
  if (!is.element(method, 
                  c("is", "nomis", "ccmn", "ruv2",
                    "ruvrand","ruvrandclust"))
  ) {
    stop("Invalid normalization method")
  }
  
  # Reference vector should be a vector
  if (!is.null(isvec)) {
    if (class(isvec) %in% c("data.frame", "list", "matrix")) {
      stop("Reference should be a vector")
    }
  }
  
  #factors should be a vector
  if (!is.null(factors)) {
    if (class(factors) %in% c("data.frame", "list", "matrix")) {
      stop("factors should be a vector")
    }
  }
  
  
  # If there is no isvec is given, get them to enter the reference vector
  if (( method == "is") & is.null(isvec)) {
    stop("Please enter the reference vector")
  }
  
  # # For the 3 ruv methods
  # # If no k, get them to enter it handled within each after plotk
  # if (substr(method,1,3) == "ruv" & is.null(k)) {
  #   stop("Please enter the number of unwanted variation factors")
  # }
  # If there is no qcmets, get them to enter it
  if ((substr(method,1,3) == "ruv" | method == "ccmn"| method == "nomis") & is.null(qcmets)) {
    stop(
      paste("Please enter a logical vector indicating",
            "quality control metabolites"
      )
    )
  }
  
  # If no factors are provided and method is ccmn, get user to enter it
  if((method == "ccmn") & is.null(factors)){
    stop(
      paste("Please enter 'factors' to be used for this method")
    )
  } else if ((method != "ccmn") & !is.null(factors))
    stop(paste("'factors' are only used in the ccmn method"))
  
  
  # if (method == "ccmn"){
  #   warning(paste("The ccmn method uses the grouping structure in the normalisation method, 
  #                 therefore, may not be used for those unsupervised 
  #                 methods where the biological factors must be treated as unknown."))      
  # }
  
  if (method == "ruv2"){
    # stop("The ruv2 method is designed for identifying 
    #         differentially abundant metabolites. Use LinearModelFit() 
    #         directly with ruv2=TRUE")     
    out<-LinearModelFit(featuredata, factormat=factormat, ruv2=TRUE,
                        k=k, qcmets=qcmets, ...)
    return(out)
  } else{
    
    # Get normalisation vector according to the method and inputs, and
    # remove groups and internal standard for data processing
    
    # Rename data for processing
    pre_norm <- featuredata
    
    if (method == "is") {
      norm_vector <- isvec
    }
    
    
    if (!is.null(qcmets)){
      ncvec<-logical(ncol(pre_norm))
      ncvec[qcmets]<-TRUE      
    }
    
    if (method == "ccmn") {
      norm_data <- t(normalize(t(pre_norm), "crmn", 
                               standards=ncvec, ncomp=ncomp,
                               factor=model.matrix(~-1+factors),
                               lg=FALSE)
      )
    } else if (method == "nomis") {
      norm_data <- t(normalize(t(pre_norm), "nomis", standards=ncvec, 
                               lg=FALSE)
      )
    } else if (method == "ruvrand") {
      norm_out<- RUVRand(Y=data.matrix(pre_norm), 
                         ctl=qcmets,lambda=lambda, k=k,plotk=plotk)
    } else if (method == "ruvrandclust") {
      norm_out_pre<- RUVRand(Y=data.matrix(pre_norm), 
                             ctl=qcmets,lambda=lambda, k=k,plotk=plotk)
      norm_out<-RuvRandIter(RUVRand=norm_out_pre, maxIter=maxIter, 
                            wUpdate=nUpdate, lambdaUpdate=lambdaUpdate, p=p)
    } else if (method == "is"){
      norm_data <- sweep(pre_norm, 1, norm_vector, "-")
    }
    
    
    # Adjust data output for different types of Nrmalisation to capture all information
    if(substr(method,1,3) == "ruv")
      outdata <- norm_out$newY
    else 
      # make dataframe for output
      outdata <- data.frame(norm_data)
    
    
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
    if(substr(method,1,3) == "ruv"){
      output <- list()
      output$featuredata <- outdata
      output$uvdata <- norm_out$UVcomp
      if (!is.null(norm_out$opt))
        output$opt<-norm_out$opt
      
    }else {
      output <- list()
      output$featuredata <- outdata
      
    }
    output$sampledata <- NULL
    output$metabolitedata <- NULL
    return(structure(output, class="alldata"))
  }
  
}


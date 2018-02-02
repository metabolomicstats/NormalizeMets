#' Generate Report
#' 
#' Generates an interactive report based on basic user input. The user can choose up to 3 normalisation methods that will
#' be compared to the unadjusted data using various diagnostics to assess the normalisation. Guidance on choosing criteria 
#' is also provided.
#' 
#' @param featuredata featuredata A data frame in the featuredata format. 
#' This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' @param sampledata A dataframe with sample information matching featuredata.
#' @param metabolitedata A dataframe with metabolite information matching featuredata.
#' @param normmeth A list of up to 3 normalisation methods. Must be one of "is",
#' "ccmn", "nomis", "ruv2", "ruvrand", "rlsc", "median", "mean", "sum". 
#' For combined methods, the list should consist of
#' vector with entries corresponding in order to the 2 methods to beused jointly.
#' @param factorOI factor of interest to be used, should correpond to column number or 
#' column name in sampledata corresponding to the factor of interest for the analysis.
#' @param gfactor A vector indicating the groups that need to be explored in the plots
#' @param covars names of the other covariates to be included when fitting the linear 
#' model for biomarker identification. Should correspond to column name in sampledata
#' @param missingvals The method to be used for removing missing values. Should be either "knn" or "replace".
#' @param logTrans A logical indication whether the data is to be log transformed.
#' @param k k Number of factors of unwanted variation to be included in the
#' "\code{ruv}" models.
#' @param fitintercept A logical indication whether an intercept component should be fitted in the linear model.
#' @param isvec A vector of internal standards to be used with the method
#' "\code{is}".
#' @param qcmets A vector indicating which metabolites should be used as the
#' internal, external standards or other quality control metabolites 
#' in the "\code{ruv}" models, or as multiple internal
#' standards in the "\code{ccmn}" and "\code{nomis}" methods. 
#' @param rlsc.sampledata For the "\code{rlsc}" method, a dataframe that contains sample specific information. 
#' Unique sample names should be provided as row names. For this function, 
#' this should have, the batch number, the class and the run order, 
#' with column names 'batch', 'class' and 'order' respectively. 
#' For the QCs samples, 'class' should be allocated as 0.
#' @param ccmn.factor For the ccmn method. A vector describing biological factors.
#' @param volcano.yrange In the volcano plot, a numeric for the maximum y value (scale of y-axis is -log(p-value)), can only
#' be set to a value as big as the maximum y-value in the plots.
#' @param scaling.refvec A reference vector for the scaling method
#' @param reportName The name that should be used to save the report.
#' @param ... Arguments to be passed onto other methods.
#' 
#' @author Alysha M De Livera, Gavriel Olshansky
#' 
#' @example 
#'      # generate an example reoprt using mixdara
#'      GenerateReport()
#'      
#' @export GenerateReport
GenerateReport <- function(featuredata = NULL, sampledata = NULL, metabolitedata = NULL,
                           normmeth = list(method1 = c("nomis"),
                                           method2 = c("ccmn"),
                                           method3 = c("ruv2")),
                           factorOI = NULL, covars = NULL, gfactor = NULL,
                           missingvals = c("knn","replace", "none"),
                           logTrans = TRUE, k = NULL, fitintercept = TRUE,
                           isvec = NULL, qcmets = NULL, rlsc.sampledata = NULL, ccmn.factor = NULL,
                           volcano.yrange = NULL,scaling.refvec = NULL,
                           reportName = "General_Report",
                           ...){
  
  # if featuredata is null, generate the example report
  if (is.null(featuredata)){
    cat("An example report will be generated using mixdata! \n")
    
    data("mixdata", envir = environment())
    featuredata <- mixdata$featuredata
    sampledata <- mixdata$sampledata
    metabolitedata <- mixdata$metabolitedata
    
    k <- 3
    missingvals = "knn"
    factorOI <- "type"
    gfactor<- "batch"
    ccmn.factor <- "batch"
    #qcmets<-which(metabolitedata$type=="IS")
    qcmets <- "type"
  }
  
  #Featuredata needs column names
  if(is.null(colnames(featuredata)))
    stop ("column names are missing from featuredata")
  

  # first make sure all input is correct..
  missingvals<-match.arg(missingvals)
  All_Norm_Methods <- list(QCmets = c("is","nomis","ccmn","ruv2", "ruvrand","ruvrandclust"),
                           QCsamples = c("rlsc"),
                           Scaling =  c("median","mean","sum","ref"))
  AllNormMethods <- c("is","nomis","ccmn","ruv2", "ruvrand","ruvrandclust","rlsc","median","mean","sum","ref")
  
  Nmethods <- length(normmeth)
  
  # make sure no more than 3 methods entered
  if (Nmethods > 3){
    stop("Enter no more than 3 normalisation methods")
  }
  
  normtype = c()
  # methods entered are recognised
  for (ii in 1:Nmethods){
    ncurrent <- length(normmeth[[ii]])
    if (ncurrent > 2){
      stop("maximum number of methods for the NormCombined is 2")
    } else if (ncurrent == 2){
      for (jj in 1:2){
        if (!(normmeth[[ii]][jj] %in% AllNormMethods))
          stop(paste("Unrecognised Normalisation method entered, choose from: is",AllNormMethods[-1], sep =", "))
        normtype[ii] <- "NormCombined"
      }
    } else if (ncurrent == 1){
      if (normmeth[[ii]][1] %in% All_Norm_Methods$QCmets){
        normtype[ii] <- "NormQcmets"
      } else if (normmeth[[ii]][1] %in% All_Norm_Methods$QCsamples){
        normtype[ii] <- "NormQcsamples"
      } else if (normmeth[[ii]][1] %in% All_Norm_Methods$Scaling){
        normtype[ii] <- "NormScaling"
      } else {
        stop(paste("Unrecognised Normalisation method entered, choose from: is",AllNormMethods[-1], sep =", "))
      }
    }
  }
  
  cat("The report is being generated, this may take a few moments.....")
  
  template_location <- system.file("rmd/General_Report_Template.Rmd", package = "NormalizeMets")

  
  # Generate the report
  rmarkdown::render(input = template_location, output_dir = getwd(), output_file = paste0(reportName,".html"), output_options = list(self_contained = TRUE),
                    params = list(featuredata = featuredata, sampledata = sampledata, metabolitedata = metabolitedata,
                                  normmeth = normmeth, normtype = normtype, Nmethods = Nmethods,
                                  factorOI = factorOI, covars = covars, gfactor = gfactor,
                                  logTrans = logTrans, missingvals = missingvals, k = k, fitintercept = fitintercept,
                                  isvec = isvec, qcmets= qcmets, rlsc.sampledata = rlsc.sampledata,
                                  ccmn.factor = ccmn.factor,
                                  volcano.range = volcano.yrange,
                                  scaling.refvec = scaling.refvec, default=TRUE))
  
  
  
  
  
  
  
}
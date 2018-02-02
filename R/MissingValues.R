#' Missing value replacement
#' 
#' Missing value imputation for metabolomics data matrices 
#' 
#' 
#' @param featuredata A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' See NormalizeMets Vignette for details.
#' @param sampledata   A dataframe with sample information matching featuredata.
#' @param metabolitedata A dataframe with metabolite information matching featuredata.
#' @param feature.cutoff A value between zero and one. Used to exclude features that have 
#' a large proportion of missing values. If the proportion of
#' missing values is equal to or more than the feature.cutoff, that feature will be deleted.
#' @param sample.cutoff A value between zero and one. Used to exclude samples that have 
#' a large proportion of missing values. If the proportion of missing
#' values is equal to or more than the sample.cutoff in any row, that whole sample
#' will be deleted.
#' @param method Missing value replacement method. Should be either "knn" (the kth nearest neighbour algorithm),
#'  "replace" (replacing by half the minimum detectable signal ), or "none".
#' @param k The number of nearest neighbours to be used in the knn algorithm
#' @param  featuremax.knn For the knn algorithm. The maximum proportion of missing data allowed in 
#' any feature. For any features with more than featuremax.knn proportion missing,
#' missing values are imputed using the overall mean per sample.
#' @param samplemax.knn For the knn algorithm. The maximum proportion of
#' missing data allowed in any sample. If any sample has more than 
#' samplemax.knn missing data, the program halts and reports an error.
#' @param seed For the knn algorithm for very large matrices. An integer, denoting state for random number generation in R.
#' @param saveoutput A logical indicating whether the output should be saved.
#' If \code{TRUE}, the results will be saved as a csv file.
#' @param outputname The name of the output file if the output has to be saved.
#' @return The output is an object of class
#' \code{\link[NormalizeMets:alldata]{alldata}}.
#' @author Alysha M De Livera, Gavriel Olshansky
#' @examples
#' 
#'     data(alldata_eg)
#'     featuredata_eg<-alldata_eg$featuredata
#'     sampledata_eg<-alldata_eg$sampledata
#'     metabolitedata_eg<-alldata_eg$metabolitedata
#'     logdata <- LogTransform(featuredata_eg)
#'     
#'     imp <-  MissingValues(logdata$featuredata,sampledata_eg,metabolitedata_eg,
#'                       feature.cutof=0.8, sample.cutoff=0.8, method="knn")
#'     imp
#'     dataview(imp$featuredata)                       
#' 
#' @export MissingValues
MissingValues <-function(featuredata,
                         sampledata=NULL,
                         metabolitedata=NULL,
                         feature.cutoff=0.8, 
                         sample.cutoff=0.8,
                         method=c("knn","replace","none"),
                         k=10,
                         featuremax.knn=0.8,
                         samplemax.knn=0.8,
                         seed=100, 
                         saveoutput=FALSE, outputname="nomissing")  #GO 10/5
{
  
  if(is.null(featuredata))
    stop("featuredata is missing")
  if(is.null(metabolitedata))
    metabolitedata<-data.frame(names=colnames(featuredata))
  if(is.null(sampledata))
    sampledata<-data.frame(names=rownames(featuredata))
  
  # sampledata should be a dataframe
  if (!class(sampledata) %in% c("data.frame")) 
    stop("sampledata should be a dataframe")
  
  # metabolitedata should be a dataframe
  if (!class(metabolitedata) %in% c("data.frame")) 
    stop("metabolitedata should be a dataframe")
  
  # sampledata should be a dataframe
  if (!class(sampledata) %in% c("data.frame")) 
    stop("sampledata should be a dataframe")
  
  method<-match.arg(method)
  
  if (method!="none"){
    if (length(method)==2 | method != "knn" & method != "replace")
      print("Method should either be knn or replace")
    set.seed(seed)
    vars <- colnames(featuredata)[1:length(colnames(featuredata))] 
    var_count <- length(vars)
    row_count <- dim(featuredata)[1]
    matrix_min <- min(featuredata, na.rm = TRUE)
    write(" -> Checking features...", "")
    
    # Remove features that have >= feature.cutoff proportion missing
    mt_cols <- which(colSums(is.na(featuredata)) >= row_count*feature.cutoff)
    if (!is.null(feature.cutoff)){
      if (length(mt_cols) > 0) {
        pretrim_data <- featuredata[, -c(mt_cols)]
        metabolitedata <- metabolitedata[-c(mt_cols),]
        write(paste("Features with ", feature.cutoff," proportion missing are removed.",sep=""),"") #GO 23/5had problem writing
        
      } else {pretrim_data <- featuredata}     #GO 16/5 
    }
    else {
      pretrim_data <- featuredata
    }
    
    write(" -> Checking samples...", "")
    # Remove samples that have > sample.cutoff proportion missing
    mt_rows <- which(rowSums(is.na(pretrim_data)) >= var_count*sample.cutoff)
    if (!is.null(sample.cutoff)){
      if (length(mt_rows) > 0) {
        pretrim_data <- pretrim_data[-c(mt_rows), ]
        sampledata <- sampledata[-c(mt_rows), ]
        write(paste("Samples with", sample.cutoff,"proportion missing are removed.",sep=" "),"") #GO 23/5 had problem writing
      }
    }
    
    pretrim_data<-as.matrix(pretrim_data)
    if (method=="knn")
    {
      featuredata<-t(impute.knn(data=t(pretrim_data) ,
                                k = k, 
                                rowmax = featuremax.knn*100, 
                                colmax = samplemax.knn*100, 
                                maxp = 1500, rng.seed=seed)$data)
      
    } else if (method == "replace"){
      if (matrix_min < 0) {
        stop("The data contains negative values.")
      }
      pretrim_data[which(is.na(pretrim_data))]<- matrix_min/2
      featuredata<-pretrim_data
      
    }
  }
  alldata <- c()
  alldata$featuredata <- featuredata
  alldata$sampledata <- sampledata
  alldata$metabolitedata <- metabolitedata
  
  if (saveoutput) {                                # GO 5/10
    write.csv(featuredata, paste(c(outputname,"_featuredata.csv"), collapse = ""))
    write.csv(sampledata, paste(c(outputname,"_sampledata.csv"),collapse=""))
    write.csv(metabolitedata, paste(c(outputname,"_metabolitedata.csv"),collapse = ""))
  }
  
  return(structure(alldata, class = "alldata"))
}


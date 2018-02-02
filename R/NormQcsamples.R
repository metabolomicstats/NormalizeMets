#' Normalisation methods based on quality control samples
#' 
#' This function is based on the quality control sample based robust LOESS 
#' (locally estimated scatterplot smoothing) signal correction (QC-RLSC) method 
#' as described by Dunn \emph{et al}. (2011) and impletemented 
#' in statTarget: Luan H (2017).
#' 
#' @param featuredata featuredata A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#'  Unique sample names should be provided as row names.
#' @param sampledata A dataframe that contains sample specific information. 
#' Unique sample names should be provided as row names. For this function, 
#' this should have, the batch number, the class and the run order, 
#' with column names 'batch', 'class' and 'order' respectively. 
#' For the QCs samples, 'class' should be allocated as 0.
#' @param method A character string indicating the required normalization
#' method. In this case, only "\code{rlsc} method as described by Dunn \emph{et al}. (2011). 
#' @param span The smooting parameter. The default is 0.75. 
# if the QCspan is set at 0, the generalised cross-validation will be performed to avoid
# overfitting the observed data.
#' @param deg The degree for the polynomial fit. Must be either 0, 1, or 2. Defaults to degree 2 for local polynomial fits.
#' @param lg A logical indicating whether the final normalized data matrix needs to be log transformed
#' after rlsc method. Defaults to TRUE.
#' @param saveoutput A logical indicating whether the normalised data matrix and related plots
#' should be saved.
#' @param outputname The name of the output file if the output has to be saved.
#' @param ... Extra input for the \code{statTarget::shiftCor} function.
#' @return The result is an object of class \code{\link[NormalizeMets:alldata]{alldata}}. 
#' @seealso \code{statTarget::shiftCor}.
#' @author Alysha M De Livera, Gavriel Olshansky
#' @references  Luan H (2017). statTarget: Statistical Analysis of Metabolite Profile. 
#' R package version 1.6.0, https://github.com/13479776/statTarget.
#' 
#' Dunn,W.B., Broadhurst,D., Begley,P., Zelena,E., Francis-McIntyre,S., Anderson,N., Brown,M., 
#' Knowles,J.D., Halsall,A., Haselden,J.N. et al. (2011) Procedures for 
#' large-scale metabolic profiling of serum and plasma using gas chromatography and 
#' liquid chromatography coupled to mass spectrometry. Nat. Protoc., 6, 1060-1083
#' 
#' @examples
#' 
#' ##Reading the data
#' data(Didata)
#' NormQcsamples(sampledata=Didata$sampledata[order(Didata$sampledata$order),],
#'               featuredata=Didata$featuredata[order(Didata$sampledata$order),],
#'               saveoutput=FALSE)
#'               
#'     
#' @export NormQcsamples
NormQcsamples <-function(featuredata,
                         sampledata,  
                         method=c("rlsc"),
                         span=0, deg=2, lg=TRUE,saveoutput=FALSE, 
                         outputname="qcsample_results", ...)
{
  # Match the method
  method <- match.arg(method)
  
  # If method is not one of the above listed, then stop
  if (!is.element(method, 
                  c("rlsc"))
  ) {
    stop("Invalid normalization method")
  }
  
  alldatacheck(featuredata=featuredata, sampledata=sampledata)
  
#   if (length(is.na(featuredata))!=0){ #shiftCor_fn let 0s equal missing
#     featuredata[which(is.na(featuredata))]<-0 
#     message("\nWarning: The zeros are assigned missing") #in the shiftCor_fn
# #    print("Warning: The zeros are assigned missing")
#   }
  
  
  featuredata_new<-data.frame(name=row.names(featuredata),featuredata,stringsAsFactors=FALSE)
  #  featuredata_new<-editcolnames(featuredata_new)
  sampledata_new<-data.frame(sample=row.names(sampledata),sampledata,stringsAsFactors=FALSE)
  
  
  #Allocate all sampledata colnames to lower
  colnames(sampledata_new)<-tolower(colnames(sampledata_new))
  
  if (!('batch' %in% colnames(sampledata_new)))
    stop(" 'batch' is missing")
  
  if (!('class' %in% colnames(sampledata_new)))
    stop(" 'class' is missing")
  
  if (!('order' %in% colnames(sampledata_new)))
    stop(" 'order' is missing")
  #Order must start from 1 and monotonically increasing
  if (sampledata_new$order[1]!=1 |  !all(sampledata_new$order == cummax(sampledata_new$order)))
      stop("Sort your data increasing order starting from 1")
  
  if (any(is.na(sampledata_new$class))==TRUE) {
    stop("There are missing values in 'class'.")
  }
  
  if (length(sampledata_new$class==0)==0) {
    stop("There are no QC samples (class==0) in your data.")
  }
  
  featuredata_new$name[which(sampledata_new$class==0)]<-paste("QC", featuredata_new$name[which(sampledata_new$class==0)],sep="_")
  sampledata_new$sample[which(sampledata_new$class==0)]<-paste("QC", sampledata_new$sample[which(sampledata_new$class==0)],sep="_")
  
  sampledata_new$class[which(sampledata_new$class==0)]<-NA
  
  #Transpose for the function
  
  #Creating .csv 
  feature_save<-t(featuredata_new)
  colnames(feature_save)<-feature_save[1,]
  feature_save<-feature_save[-1,]
  write.csv(feature_save, "temp_feature.csv")
  #  write.table(t(featuredata_new), "temp_feature.csv",col.names=FALSE)
  write.csv(sampledata_new, "temp_sample.csv", row.names = FALSE)
  samPeno <- paste(getwd(),"temp_sample.csv", sep="/")
  samFile <- paste(getwd(),"temp_feature.csv", sep="/")
  
  outdata<-shiftCor_fn(samPeno,
                           samFile,
                           QCspan = span, 
                           degree = deg,
                           saveoutput=saveoutput, 
                           outputname=outputname,
                           lg=lg,
                           ...)
  
  file.remove("temp_feature.csv")
  file.remove("temp_sample.csv")
  
#  outdata<-list()
  # outdata$sampledata<-data.frame(sampledata,outdata_pre[,c(1:2)])
  # outdata$metabolitedata<-NULL
  # outdata$featuredata<-as.matrix(outdata_pre[,c(3:dim(outdata_pre)[2])])
  # rownames(outdata$featuredata)<-rownames(outdata$sampledata)
  # outdata$featuredata<-editcolnames(outdata$featuredata)
  
# Save data if required
# Generate output name in a csv format 
  if (saveoutput) {   #GO - 20/17
    write.csv(outdata,
              if (!is.null(outputname)) {
                paste(c(outputname, ".csv"), collapse="")
              } else {
                paste(c("normalized_", method, ".csv"), collapse="")
              }
    )
  }

  
  return(structure(outdata, class="alldata"))
  
  
}


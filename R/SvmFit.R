#'  support vector machine
#' 
#' Classification using support vector machine (svm) algorithm
#'
#' @param featuredata featuredata A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' See NormalizeMets Vignette for details.
#' @param groupdata A vector with group names.
#' @param cost of constraint violation, defaults to 1.
#' @param gamma parameter used for the kernel
#' @param kernel The kernel used. The default is the radial basis function with type C-classification.
#' @param crossvalid A logical indicating whether cross-validation needs to be conducted
#' @param k An integer specifying the k-fold cross-validation. Default is set to 5.
#' @param tune A logical with the default set to FALSE. If TRUE, a grid search will be 
#' conducted to tune the hyperparameters, over parameter ranges supplied by the user.
#' @param pred whether the predictions should be made
#' @param pfeaturedata The test dataset for the predictions. The default is featuredata
#' @param pgroupdata The test groupdata for the predictions. The default is groupdata
#' @param rocplot A logical indicating whether a receiver operating characteristic curve needs to be plotted, 
#' along with the area under the curve (AUC) printed.
#' @param saveoutput A logical indicating whether the outputs should be saved in the format
#' \code{\link[e1071]{write.svm}}
#' @param main Plot title.
#' @param outputname The name of the output file if the output has to be saved.
#' @param ... Arguments to be passed on to other methods.
#' @return If tune=FALSE, an object of class "svm" \code{\link[e1071]{svm}} 
#' containing the fitted model or if tune=TRUE, 
#' an object of class \code{\link[e1071]{tune}} 
#' @author Alysha M De Livera, Gavriel Olshansky
#' @seealso \code{\link[e1071]{tune}}
#' @examples
#'  data(alldataC)
#' SvmFit(featuredata=alldataC$featuredataC, 
#'        groupdata=alldataC$groupdataC,
#'        crossvalid=TRUE,
#'        k=5,
#'        rocplot = TRUE)
#' 
#' @export SvmFit
SvmFit <- function(featuredata, groupdata, 
                   kernel= "radial",
                   cost=1,
                   gamma=NULL,
                   crossvalid=TRUE,
                   k=5,
                   tune=FALSE,
                   pred=TRUE,
                   pfeaturedata=featuredata,
                   pgroupdata=groupdata,
                   rocplot=TRUE,
                   saveoutput=FALSE,
                   outputname="svm",
                   main=NULL,
                   ...){
  
  
  # Enter groupdata as a vector
  if (is.null(groupdata))
    stop("Enter groupdata")
  else if (class(groupdata) %in% c("data.frame", "list", "matrix")) {
    stop("groupdata should be a vector")
  }
  
  groupdata<-as.factor(groupdata)
  inputdata<-data.frame(group=groupdata,featuredata)
  
  if (crossvalid)
    cross=k
  else
    cross=0
  
  svmfit<-svm(group ~ ., data = inputdata, cross=cross,
              type= "C-classification",
              kernel=kernel,
              ...)
  
  if (tune)
    tuned<-tune.svm(group ~ ., data = inputdata,
                    gamma = gamma, cost = cost, 
                    ...)
  
  
  if(pred){
    if(tune)
      stop("predictions cannot be made when tune is set to TRUE")
    else {
      preds<-predict(svmfit,pfeaturedata,type="class")
      if (rocplot)
        plot(roc(predictions=preds, labels=pgroupdata),
             col="blue",
             lwd=2,main=main)
      text(0.5,0.5, paste("AUC=", round(auc(roc(predictions=preds, labels=pgroupdata)),2)),
           col="brown")
      
    }
    
  }
  
  
  if (tune)
    return(tuned)
  else
    return(svmfit)

  if (saveoutput)
    write.svm(svmfit, svm.file = paste(outputname,"classifier.svm"), 
              scale.file = paste(outputname,"classifier.scale"),
              yscale.file =  paste(outputname,"classifier.yscale"))

  
  
}


  
#' p-value Histogram
#' 
#' Make a p-value Histogram of results
#' 
#' @param lpvals A list with vectors of p-values
#' @param normmeth A vector with the normalization method used corresponding in order to the 
#' data supplied to be displayed on the plot.
#' @param saveplot A logical indication whether to save the plot produced.
#' @param plotname Name of the output file if the file is to be saved. This is
#' the general name for all the graphs and the specific type prefix will be
#' added automatically.
#' @param savetype The required format for the plot to be saved in. Threre is a
#' choice of \code{"png","bmp","jpeg","tiff","pdf"} type files.
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ylim y-axis limit
#' @param xlim x-axis limit
#' @param col	a colour to be used to fill the bars. The default of NULL yields unfilled bars.
#' @param ...	Other parameters for the \code{\link[graphics]{hist}} function.
#' @author Alysha M De Livera, Gavriel Olshansky
#' @examples 
#' data("alldata_eg")
#' featuredata_eg<-alldata_eg$featuredata
#' dataview(featuredata_eg)
#' sampledata_eg<-alldata_eg$sampledata
#' dataview(sampledata_eg)
#' metabolitedata_eg<-alldata_eg$metabolitedata
#' dataview(metabolitedata_eg)
#'
#' logdata <- LogTransform(featuredata_eg)
#' dataview(logdata$featuredata)
#' imp <-  MissingValues(logdata$featuredata,sampledata_eg,metabolitedata_eg,
#'                      feature.cutof=0.8, sample.cutoff=0.8, method="knn")
#' dataview(imp$featuredata)
#'
#' #Linear model fit using unadjusted data
#' factormat<-model.matrix(~gender +Age +bmi, sampledata_eg)
#' unadjustedFit<-LinearModelFit(featuredata=imp$featuredata,
#'                              factormat=factormat,
#'                              ruv2=FALSE)
#' unadjustedFit
#' 
#' #Linear model fit using `is' normalized data 
#' Norm_is <-NormQcmets(imp$featuredata, method = "is", 
#'                     isvec = imp$featuredata[,which(metabolitedata_eg$IS ==1)[1]])
#' isFit<-LinearModelFit(featuredata=Norm_is$featuredata,
#'                      factormat=factormat,
#'                      ruv2=FALSE)
#' isFit
#'
#' #Linear model fit with ruv-2 normalization
#' ruv2Fit<-LinearModelFit(featuredata=imp$featuredata,
#'                        factormat=factormat,
#'                        ruv2=TRUE,k=2,
#'                        qcmets = which(metabolitedata_eg$IS ==1))
#' ruv2Fit
#'
#' #Exploring metabolites associated with age
#' lpvals_age<-list(unadjusted=unadjustedFit$p.value[,"Age"],
#'                 is=isFit$p.value[,"Age"],
#'                 ruv2=ruv2Fit$p.value[,"Age"])
#'
#' ComparePvalHist(lpvals = lpvals_age,ylim=c(0,40), 
#'      normmeth = c("unadjusted","is","ruv2"))
#' 
#' @export ComparePvalHist
ComparePvalHist <- function(lpvals=NULL,
                            normmeth=NULL, 
                            saveplot = FALSE,
                            savetype= c("png","bmp","jpeg","tiff","pdf"),
                            xlab = "P-Values",
                            ylab= "Frequency",
                            ylim=NULL,
                            xlim=c(0,1),
                            col="grey",
                            plotname = "PvalHistComp",
                            ...){
  if (is.null(lpvals))
    stop("Enter a list of p-values")
  nplots <- length(lpvals)      # get total number of plots to compare
  
  #get format for saving the file if needed
  savetype <- match.arg(savetype)
  
  #save plot if needed                             #GO 18/7
  
  if (saveplot == TRUE) { 
    savef <- match.fun(savetype)
    if (savetype == "pdf") {
      savef(paste(c(plotname,".",savetype),collapse = ""),width = 10*(floor(nplots/4)+1), height = 6*(nplots-floor(nplots/4)*2))
    }else {
      savef(paste(c(plotname,".",savetype),collapse = ""),width = 1000*(floor(nplots/4)+1), height = 600*(nplots-floor(nplots/4)*2))
    }
  }
  
  if (nplots <= 2)
    op <- par(mfrow=c(nplots,1))
  else if (nplots ==3)
    op <- par(mfrow=c(3,1))
  else if (nplots <=4)
    op <- par(mfrow=c(2,2))
  else
    stop ("A maximum of four plots are compared.")
  
  if(is.null(normmeth) | length(normmeth) != nplots){
    for (ii in 1:nplots){
      normmeth <- c(normmeth,paste("Method",ii))
    }
  }
  

    
  for (jj in 1:nplots){
    hist(lpvals[[jj]], 
         main = normmeth[[jj]], 
         xlab = xlab ,
         ylab =ylab,
         xlim=xlim,
         ylim=ylim,
         col=col,
         ...)
  }
  
  if (saveplot == TRUE) {
    dev.off()
  }
  
  par(op)
}

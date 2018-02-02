#' Compare Volcano plots
#' 
#' Produces a volcano plot that can be sued to compare between different normalisation methods.
#' 
#' @param lcoef A list of vectors of coefficients with metabolite names, each vector corresponding to a 
#' different noramlization method.
#' @param lpvals A list of vector of corresponding p-values.
#' @param normmeth A vector with the normalization method used corresponding in order to the 
#' data supplied to be displayed on the plot.
#' @param plimit A numeric indicating the p value cutoff. The default is set to
#' 0.05.
#' @param coeflimit A numeric indicating the lower fold cutoff. The default is
#' set to 2.
#' @param yrange A numeric for the maximum y value (scale of y-axis is -log(p-value)), can only
#' be set to a value as big as the maximum y-value in the plots.
#' @param negcontrol A vector with the names of the metabolites used as negative controls,
#' to be coloured differently.
#' @param poscontrol A vector with the names of the metabolites used as positive controls,
#' to be coloured differently.
#' @param xlab \emph{x}-axis label.
#' @param ylab \emph{y}-axis label.
#' @param labelunderlim A logical indicating whether to label points that are not significant.
#' @param labelsig A logical indicating whether all significant points should be labeled.
#' @param saveinteractiveplot A logical indication whether the interactive plot produced should
#' be saved as a \code{.html} file.
#' @param interactiveplotname A character string indicating the name to be used for saving the
#' interactive plot.
#' @param ... Arguments to VolcanoPlot function
#' @seealso \code{\link{VolcanoPlot}}
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
#' lcoef_age<-list(unadjusted=unadjustedFit$coefficients[,"Age"],
#'                is_age=isFit$coefficients[,"Age"],
#'                ruv2_age=ruv2Fit$coefficients[,"Age"])
#' lpvals_age<-list(unadjusted=unadjustedFit$p.value[,"Age"],
#'                 is=isFit$p.value[,"Age"],
#'                 ruv2=ruv2Fit$p.value[,"Age"])
#'
#' negcontrols<-metabolitedata_eg$names[which(metabolitedata_eg$IS==1)]                   
#'
#' CompareVolcanoPlots(lcoef=lcoef_age, 
#'                    lpvals_age, 
#'                    xlab="Coef",
#'                    negcontrol=negcontrols)
#'
#'@export CompareVolcanoPlots
CompareVolcanoPlots <- function(lcoef, lpvals, normmeth = NULL, 
                                plimit=0.05, 
                                coeflimit=1,
                                yrange=NULL,
                                negcontrol = NULL, poscontrol=NULL, 
                                xlab='Coefficients',
                                ylab='-log(p-value)', labelunderlim = FALSE, 
                                labelsig=FALSE,
                                saveinteractiveplot = FALSE,
                                interactiveplotname="interactiveVolcanPlot",
                                ...){
  
  
  # make sure input dimensions are right
  if (length(lcoef) != length(lpvals)){
    stop("Dimensions of data entered must match - make sure list of folds and pvals are of
         the same length")
  }
  
  nplots <- length(lcoef)      # get total number of plots to compare
  
  # set legend parameter to pas to volcanoplot function
  currentlegend <- rep(TRUE,6)
 
  
############################################################################################################
  #Generate names to distinguish between plots - incase they are not given
  if(is.null(normmeth) | length(normmeth) != nplots){
    for (ii in 1:nplots){
      normmeth <- c(normmeth,paste("Method",ii))
    }
  }
  
  Iplots <- list(rep(plot_ly(type = 'box'),nplots))
  subp <- "subplot("             # store the subplot command as string for later use (to be able to add the 
  # required number of plots as need to add them all at once)
  
  
  #generate the seperate plots and add them to the subp string
  for (jj in 1:nplots){
    CplotOut <- VolcanoPlot(lcoef[[jj]],
                lpvals[[jj]],
                plimit=plimit,
                coeflimit=coeflimit,
                xlab=xlab,
                ylab=paste(ylab, normmeth[[jj]],sep= " ") ,
                negcontrol = negcontrol,
                poscontrol = poscontrol,
                labelunderlim= labelunderlim,
                labelsig = labelsig,
                saveinteractiveplot=FALSE,#because we don't want to save this one
                interactiveplot = TRUE,
                interactiveonly = TRUE,
                chooselegend = currentlegend,
                ...)
    
    Iplots[[jj]] <- CplotOut[[1]]
    currentlegend <- CplotOut[[2]]
        
    #set subplot title
    if (!is.null(yrange)){
      Iplots[[jj]] <- add_trace(Iplots[[jj]],x=0,y=max(yrange),mode="markers",type = "scatter",
                                marker = list(size=1, color = "black"),hoverinfo = "none",
                                showlegend = FALSE)
      Iplots[[jj]] <- layout(Iplots[[jj]])
    }
    subp<- paste(subp,"Iplots[[",jj,"]],",sep="")
    
  }
  
  subp <- paste(subp,"shareX = TRUE, shareY = TRUE, titleY = TRUE, nrows=",nplots,")",sep = "")
  
  p <- eval(parse(text = subp))
  
  p <- layout(p)
  
  
  if (saveinteractiveplot){                                               
    htmlwidgets::saveWidget(p, paste(interactiveplotname,".html",sep=""))
  }
  
  
  return(p)
  
##########################################################################################################
  
}
#' CompareRlaPlots
#' 
#' Produces within group and across group relative log abundance plots to
#' visually compare between different normalization methods
#' 
#' @param lfeaturedata A list containing data frames in the featuredata format.
#' @param groupdata A vector containing group information. 
#' @param normmeth A vector with the normalization method used corresponding in order to the 
#' data supplied to be displayed on the plot.
#' @param type A character string indicating whether within group ("\code{wg}")
#' or across group ("\code{ag}") RLA plots need to be plotted.
#' @param yrange A vector with the first entry corresponding to the minimum y-axis value and the second
#' to the maximum y-axis value to show as default on all the plots. This can be zoomed out 
#' as the plot is interactive.
#' @param plottitle The title to be displayed on the plot.
#' @param saveinteractiveplot A boolean indicating whether the interactive plto should be save as a 
#' \code{.html} file.
#' @param savenoninteractive A boolean indicating whether a \code{.png} version of the plto should be save as a.
#' @param interactivesavename A character string to be used as the filename for saving the interactive plot.
#' @param ... Other arguments to \code{\link[NormalizeMets]{RlaPlots}} function.
#' 
#' @examples 
#' data(UVdata)
#' # Not RUN due to user input; we set k=1 each and saved normalized data as uv_ruvrandclust
#' # uv_ruvrand_norm<-NormQcmets(featuredata=UVdata$featuredata,
#' #                            method="ruvrandclust",
#' #                            qcmets=which(UVdata$metabolitedata$neg_control==1),
#' #                            k=1)
#' data("uv_ruvrandclust")
#' lfeaturedata<-list(unadj=UVdata$featuredata,ruv=uv_ruvrandclust$featuredata,
#'    ruvuv=uv_ruvrandclust$uvdata)
#' #CompareRlaPlots(lfeaturedata,
#' #                groupdata=interaction(UVdata$sampledata$temperature,UVdata$sampledata$instrument),
#' #                normmeth=c("Unadjusted:", "RUVrandclust normalized:", 
#' #                           "RUVrandclust: removed uv:"),
#' #               yrange=c(-3,3))
#' 
#' @export CompareRlaPlots
CompareRlaPlots <- function(lfeaturedata, groupdata, normmeth=NULL, type=c("ag", "wg"),yrange=NULL,
                            plottitle = "RLA plots Comparison", saveinteractiveplot = FALSE,
                            savenoninteractive = FALSE,
                            interactivesavename = "RlaPlotsComp",...){
  
  
  # #make sure input dimensions are right
  # if(length(lfeaturedata) != length(groupdata)){
  #   stop("The number of groupings (length of groupdata) needs to be the same number 
  #        of datasets (length of lfeaturedata)")
  # }
  if (class(groupdata) %in% c("data.frame", "list", "matrix")) 
    stop("groupdata should be a vector")
    
  
  nplots <- length(lfeaturedata)      # get total number of plots to compare
  
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
  needlegend <- TRUE         # a legend needs to be made for the first plot
  
  for (jj in 1:nplots){
    Iplots[[jj]] <- RlaPlots(lfeaturedata[[jj]],groupdata,type = type, interactiveplot=TRUE,
                             interactiveonly=TRUE,showlegend = needlegend,...)
    #set subplot title
    Iplots[[jj]] <- layout(Iplots[[jj]],yaxis = list(title = paste(normmeth[jj], "RLA" ,sep = " ")))
    
    subp<- paste(subp,"Iplots[[",jj,"]],",sep="")
    
    needlegend <- FALSE      #legend only needed for the first plot
    
  }
  
  # set subplots range 
  for (jj in 1:nplots){
    if (!is.null(yrange)){
      Iplots[[jj]] <- layout(Iplots[[jj]],yaxis = list(range = yrange))
    }
    ## write something to set range to smallest range ....
  }
  
  subp <- paste(subp,"shareX = TRUE, titleY = TRUE, nrows=",nplots,")",sep = "")
  
  p <- eval(parse(text = subp))
  
  p <- layout(p,title = plottitle)
  
  
  if (saveinteractiveplot){                                               
    htmlwidgets::saveWidget(p, paste(interactivesavename,".html",sep=""))
  }
  #if (savenoninteractive){
  #  plotly_IMAGE(p, format = "png", out_file = paste(interactivesavename,".png",sep=""))
  #}
  
  
  return(p)

  
}


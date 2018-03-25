#' RLA plots
#' 
#' Produces within group and across group relative log abundance plots to
#' visualise a metabolomics data matrix. See De Livera \emph{et al}. 2012, 2013, and 2015 for details.
#' 
#' @param featuredata featuredata A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' See NormalizeMets Vignette for details.
#' @param groupdata A data frame or a vector with group names.
#' @param minoutlier A number indicating which samples names to show, samples
#' names will only be shown for samples with a median reading greater than the
#' number entered.
#' @param type A character string indicating whether within group ("\code{wg}")
#' or across group ("\code{ag}") RLA plots need to be plotted.
#' @param saveplot A logical indication whether to save the plot produced.
#' @param plotname Name of the output file if the file is to be saved.
#' @param savetype The required format for the plot to be saved in. Threre is a
#' choice of \code{"png","bmp","jpeg","tiff","pdf"} type files.
#' @param interactiveplot A boolean indicator whether to make an interactive plot.
#' @param interactiveonly A boolean indicating whether only an interactive plot should be returned.
#' @param saveinteractiveplot A boolean indicating whether to save the interactive plot as an html file.
#' @param interactivesavename A character string to be used as the filename for saving the interactive plot.
#' @param cols A character string with colours to be used for the box plots.
#' @param cex.axis The magnification to be used for \emph{x}- and
#' \emph{y}-labels relative to the current setting of cex.
#' @param las A numeric in 0, 1, 2, 3 denoting the style of axis labels.  See
#' \code{\link[graphics]{par}}.
#' @param keeporder A logical indicator whether to keep the original sample order or group them by groupdata.
#' @param ylim A vector containing \emph{y}-axis limits.
#' @param oma A vector giving the size of the outer margins.
#' @param xlabel Label for the x-axis
#' @param showlegend A logical indicator whether to display a legend for the plot.
#' @param ... Other graphical parameters. See \code{\link[graphics]{par}}.
#' @author Alysha M De Livera, Gavriel Olshansky
#' @references De Livera, A. M., Dias, D. A., De Souza, D., Rupasinghe, T.,
#' Pyke, J., Tull, D., Roessner, U., McConville, M., Speed, T. P. (2012a)
#' Normalising and integrating metabolomics data. \emph{Analytical Chemistry}
#' 84(24): 10768-10776.
#' 
#' De Livera, A. M., Olshansky M., and Speed, T.P. 2013. Statistical 
#' Analysis of Metabolomics Data. Methods in Molecular Biology (Clifton, N.J.) 
#' 1055 (Jan): 291-307.
#' 
#' De Livera, A. M,. Aho-Sysi M., Laurent J., Gagnon-Bartch J., 
#' Sandra, C, Simpson, J.A., and Speed, T. P. 2015. Statistical Methods 
#' for Handling Unwanted Variation in Metabolomics Data. Analytical Chemistry 87 (7). 
#' American Chemical Society: 3606-3615. 
#' @examples
#' 
#'     data(mixdata)
#'     
#'     RlaPlots(mixdata$featuredata, mixdata$sampledata[,1], 
#'     ylim = c(-2, 2), cols = c("green","purple"),cex.axis = 0.8)
#'     
#' @export RlaPlots
RlaPlots <- function(featuredata, groupdata, minoutlier = 0.5, type=c("ag", "wg"), 
                     saveplot=FALSE, plotname = "RLAPlot", savetype= c("png","bmp","jpeg", "tiff","pdf"),
                     interactiveplot=TRUE, interactiveonly = TRUE, saveinteractiveplot = FALSE,
                     interactivesavename = "interactiveRlaPlot", cols=NULL, cex.axis=0.7, las=2, 
                     keeporder = FALSE, ylim=NULL, oma=c(3, 3, 3, 5) + 0.1, 
                     xlabel="Samples", showlegend = TRUE,...)
{
  
  

  type <- match.arg(type)
  savetype <- match.arg(savetype)
  
  # make sure groupdata is in a data frame fromat
  groupdata <- as.data.frame(groupdata)
  
  # get unique groups
  groups <- factor(groupdata[, 1], levels = unique(groupdata[, 1]))
  unique.groups <- levels(groups)
  
  # Adjust logical
  if (interactiveplot == FALSE){
  	interactiveonly <- FALSE
  }

 
  # Get the median and standardise the data matrix, reordering to cluster by group
  # also re-order group data to match new data matrix (out_data)
  
  groupdata[,1] <- as.character(groupdata[,1])
  
  # Prepare to store new ordered data
  out_data<-data.frame()
  out_group <- c()
  group_len <- c()     #not used
  
  # Within groups
  if(type == "wg") {
    for (grp in unique.groups) {
      submat <- featuredata[which(groupdata[, 1] == grp), ]
      subgrp <- groupdata[which(groupdata[, 1] == grp), 1]
      med_vals <- apply(submat, 2, function(x) median(x,na.rm = TRUE))
      swept_mat <- sweep(submat, 2, med_vals, "-")
      out_group <- c(out_group, subgrp)
      out_data <- rbind(out_data, swept_mat)       #bind and order samples by groups
    }
    # Across groups (i.e. type == "ag")
  } else  {
    med_vals <- apply(featuredata, 2, function(x) median(x,na.rm = TRUE))
    unordered_data <- sweep(featuredata, 2, med_vals, "-")
    for (grp in unique.groups){
      subdata <- unordered_data[which(groupdata[, 1] == grp), ]
      subgrp <- groupdata[which(groupdata[, 1] == grp), 1]
      group_len <- c(group_len,length(subgrp))   #not used
      out_group <- c(out_group, subgrp)
      out_data <- rbind(out_data, subdata)       #bind and order samples by groups
    }
  }
  
  out_group <- as.data.frame(out_group)
  
  
  # get groups and dedicate group colours
  if (!keeporder){
    groups <- factor(out_group[, 1], levels = unique(out_group[, 1]))
  } else {
    out_data <- unordered_data
    out_group[,1] <- groups
  }
  
  unique.groups <- levels(groups)
  if (is.null(cols)) 
    cols <- ColList(length(unique.groups))
  box_cols <- c(rep(NA, length(rownames(out_group))))
  for (ii in 1:length(out_group[, 1])) 
    box_cols[ii] <- cols[which(unique.groups == out_group[, 1][ii])]
  
  
  #only make normal plot if required
  if(interactiveonly == FALSE){
    
    # Record Par setting and restore them when exiting
    .parold <- par(no.readonly = TRUE)
    on.exit(par(.parold))
    
  
    if (saveplot == TRUE) {
      savef <- match.fun(savetype)
      if (savetype == "pdf") {
        savef(paste(c(plotname,".",savetype),collapse = ""),width = 11, height = 9)
      }else {
        savef(paste(c(plotname,".",savetype),collapse = ""),width = 960, height = 720)
      }
    }
  
    #Find outlying samples and prepare to plot their names..
    med_vals <- apply(out_data, 1, function(x) median(x,na.rm = TRUE))
    sample_n <- vector(,length(med_vals))
    for (i in 1:length(med_vals)){
      temp.i <- i/10
      if (temp.i - round(temp.i) == 0){
        sample_n[i] <- rownames(out_data)[i]
      } else {
        sample_n[i] <- " "
      }
    }
    
    par(mar = oma)
  
    boxplot(t(out_data),
          cex.lab = 1,
          cex.axis = cex.axis,
          cex.main = 1.2,
    	    main = plotname,
          las=las,                           # label orientation
          col=box_cols,                      # colours
          ylim=ylim,                         # y-axis range                         # outer margin size
          boxlwd=0.4,
          whisklty=1,
          whisklwd=0.5,
          whiskcol="grey",
          outpch=20,
          outcex=0.1,
          show.names=FALSE,
          xlab="",
          ...
    )
    
    axis(1, 1:length(med_vals), labels=sample_n, las=2, line=(-0.5), tick=0, 
        cex.axis=cex.axis, cex.lab = 1,)
    
    title(xlab = xlabel, line = 1.8, cex.lab = 1)
    abline(h=0)
  
    par(xpd=TRUE)        #allow the legend to be drown outside the original plot
  
    # legend to be ploted differently if the plot is saved or ploted within r
    if (saveplot == TRUE) {
    
      legend("right",inset = c(1,-0.06),bg = "white", legend = unique.groups,
         col=cols, lty = c(1,1),lwd = c(2.5,2.5),cex = 1.2, horiz = TRUE)
      dev.off()
  
      } else {
      legend("left",inset = c(1,0), legend = unique.groups,
            col=cols, lty = c(1,1),lwd = c(2.5,2.5),cex = 0.7, horiz = FALSE)
      }
  }
  
  
  
  if (interactiveplot == TRUE){
    
    if (interactiveonly == FALSE){
      cat("Hit <Return> to see interactive plot: ")
      line <- readline()
    }
    
    
    total_samples <- nrow(out_data)
    boxtrace <- rep(list(),total_samples)
    
    for (i in 1:total_samples){
      boxtrace[[i]] <- list(
        y = as.numeric(out_data[i,]),
        name = as.character(rownames(out_data)[i]),
        boxpoints = 'outliers',
        marker = list(color = box_cols[i],size = 1),
        line = list(color = box_cols[i],width=1)
      )
    }
    
    p<-plot_ly(type = 'box')
    
    for (j in 1:nrow(out_data)){
      ## custom hover info not supported in plotly for boxplots at the moment
     # hover_info <- paste("Sample ref: ",rownames(out_data)[j],
     #      "</br> Median: ",median(boxtrace[[j]]$y),
     #      "</br> Variance: ",var(boxtrace[[j]]$y ),
     #      sep = "")
      
      
      p <- add_trace(p, type = 'box', y = boxtrace[[j]]$y,x=j, 
                     name = rownames(out_data)[j], #name = paste("sample: ",rownames(out_data)[j]),
                       boxpoints = boxtrace[[j]]$boxpoints, jitter = 0, whiskerwidth = 1,
                       marker= boxtrace[[j]]$marker, line = boxtrace[[j]]$line, legendgroup = box_cols[j],
                       showlegend= FALSE)
      
    }
    
    #add Legend
    if (showlegend){
      ii <- 1      #set reference for colour to use
      for (grp in unique.groups){
        p <- add_trace(p,y = 0,x=0, name = grp, hoverinfo = "none",
                      boxpoints = 'outliers', opacity = 0.5,
                      marker = list(color = cols[ii]),
                      line = list(color = cols[ii]), boxmean = FALSE,
                      legendgroup = cols[ii],showlegend = TRUE)
        ii <- ii+1
      }
    }
    
    p <- layout(p,title = plotname, xaxis = list(title = xlabel, range = c(0.5,total_samples+0.5)), yaxis = list(range = ylim))
    
    if (saveinteractiveplot){                                               
      htmlwidgets::saveWidget(p, paste(interactivesavename,".html",sep=""))
    }
      
    return(p)  
    
    
  }
  
  
}

#' Venn Diagram
#' 
#' Produces a Venn diagram showing the number of common metabolites.
#' 
#' 
#' @param lnames A list of up to three vectors, e.g. metabolite names.
#' @param group.labels A vector of reference values to be plotted, such as an
#' internal standard or sample weights.
#' @param saveplot A logical indication whether to save the plot produced.
#' @param plotname Name of the output file if the file is to be saved. This is
#' the general name for all the graphs and the specific type prefix will be
#' added automatically.
#' @param savetype The required format for the plot to be saved in. Threre is a
#' choice of \code{"png","bmp","jpeg","tiff","pdf"} type files.
#' @param main A title for the plot.
#' @param cexval The font size of the text labels.
#' @param asp The aspect ratio of the plot. A value of 1 produces a square plot
#' region.
#' @param ... Other graphical parameters. See \code{\link[graphics]{par}}.
#' @author Alysha M De Livera
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
#' lnames<- list(names(ruv2Fit$coef[,"Age"])[which(ruv2Fit$p.value[,"Age"]<0.05)],
#'              names(unadjustedFit$coef[,"Age"])[which(unadjustedFit$p.value[,"Age"]<0.05)],
#'              names(isFit$coef[,"Age"])[which(isFit$p.value[,"Age"]<0.05)])
#'
#' VennPlot(lnames, group.labels=c("ruv2","unadjusted","is"))
#'
#' @export VennPlot
VennPlot <- function (lnames, 
                      group.labels=c("A", "B", "C"),
                      saveplot = FALSE,
                      savetype= c("png","bmp","jpeg","tiff","pdf"),
                      plotname = "VennPlot",
    main="Venn Diagram", cexval=1, asp=1, ... )
{
  
  
  #get format for saving the file if needed
  savetype <- match.arg(savetype)
  
  #save plot if needed                             #GO 18/7
  
  if (saveplot == TRUE) { 
    savef <- match.fun(savetype)
    if (savetype == "pdf") {
      savef(paste(c(plotname,".",savetype),collapse = ""),width = 10, height = 10)
    }else {
      savef(paste(c(plotname,".",savetype),collapse = ""),width = 1000, height = 1000)
    }
  }
  
  
  
    # Define a list of colours to use
    cols <- ColList(7)

    # Count number of groups
    n_groups <- length(lnames)
    if (!(n_groups == 2 | n_groups == 3)) {
        stop("Only able to create Venn diagrams for 2 or 3 sets.")
    }

    # Get labels for each group
    group_labels <- group.labels[1:n_groups]
    for (ii in 1:n_groups) {
        if (length(names(lnames)) == 0) {
            group_labels[ii] <- paste(
                group_labels[ii],
                "\n(",
                length(lnames[[ii]]),
                ")",
                sep=""
            )
        } else {
            group_labels <- paste(
                names(lnames),
                "\n(",
                length(lnames[[ii]]),
                ")",
                sep=""
            )
        }
    }

    # Get values for unions
    if (n_groups == 2) {
        union_12 <- length(which(lnames[[1]] %in% lnames[[2]]))
        # Get values for single sets
        single_1 <- length(lnames[[1]]) - union_12
        single_2 <- length(lnames[[2]]) - union_12
    } else {
        union_123 <- length(
            which(
                lnames[[1]][
                    which(lnames[[1]] %in% lnames[[2]])
                ] %in% lnames[[3]]
            )
        )
        union_12 <- length(
            which(lnames[[1]] %in% lnames[[2]])
        ) - union_123
        
        union_13 <- length(
            which(lnames[[1]] %in% lnames[[3]])
        ) - union_123
        
        union_23 <- length(
            which(lnames[[2]] %in% lnames[[3]])
        ) - union_123
        
        # Get values for single sets
        single_1 <- length(lnames[[1]]) - (union_12 + union_13 + union_123)
        single_2 <- length(lnames[[2]]) - (union_12 + union_23 + union_123)
        single_3 <- length(lnames[[3]]) - (union_13 + union_23 + union_123)

    }
    
    
    
    # Produce an empty plot from 0 to 1 on both axes
    plot(c(0, 1), c(0, 1), type="n", main=main, 
        xlab="", ylab="",
        asp=asp,
        axes=FALSE
    )

    if (n_groups==2) {
    # Collect the values that will be put into the plot
        # Draw two overlapping circles
        circle_fn(0.3, 0.5, 0.25, nv=1e5, lwd=2, border=cols[1])
        circle_fn(0.6, 0.5, 0.25, nv=1e5, lwd=2, border=cols[2])
        # Add in the values to the respective sections
        text(0.25, 0.5, single_1, col=cols[1])
        text(0.45, 0.5, union_12, col=cols[4])
        text(0.65, 0.5, single_2, col=cols[2])
        # Add in group labels
        text(0.3, 0, group_labels[1], col=cols[1])
        text(0.6, 0, group_labels[2], col=cols[2])
    } else {
        # Draw three overlapping circles
        circle_fn(0.38, 0.375, 0.25, nv=1e5, lwd=2, border=cols[1])
        circle_fn(0.62, 0.375, 0.25, nv=1e5, lwd=2, border=cols[2])
        circle_fn(0.5, 0.62, 0.25, nv=1e5, lwd=2, border=cols[3])
        # Add in the values to the respective sections
        text(0.29, 0.32, single_1, col=cols[1])
        text(0.72, 0.32, single_2, col=cols[2])
        text(0.5, 0.7, single_3, col=cols[3])
        text(0.5, 0.28, union_12, col=cols[4])
        text(0.36, 0.53, union_13, col=cols[5])
        text(0.64, 0.54, union_23, col=cols[6])
        text(0.5, 0.47, union_123, col=cols[7])

        # Add in group labels
        text(0.29, 0.05, group_labels[1], col=cols[1])
        text(0.72, 0.05, group_labels[2], col=cols[2])
        text(0.5, 0.95, group_labels[3], col=cols[3])
    }

    if (saveplot == TRUE) {
      dev.off()
    }
    
}

#' Volcano plot
#' 
#' Produces a volcano plot given fold changes and p-values.
#' 
#' 
#' @param coef A vector of coefficients with metabolite names.
#' @param pvals A vector of corresponding p-values.
#' @param cexcutoff Font size of the cut-off labels.
#' @param cexlab Font size of the variable labels.
#' @param pointsize A numeric indicating the size of the points on the graph.
#' @param plimit A numeric indicating the p value cutoff. The default is set to
#' 0.05.
#' @param negcontrol A vector with the names of the metabolites used as negative controls,
#' to be coloured differently.
#' @param poscontrol A vector with the names of the metabolites used as positive controls,
#' to be coloured differently.
#' @param saveplot A logical indicator whether to save the produced plot.
#' @param plotname A character string indicating the name to be used for saving the plot.
#' @param savetype The required format for the plot to be saved in. Threre is a
#' choice of \code{"png","bmp","jpeg","tiff","pdf"} type files.
#' @param coeflimit A numeric indicating the lower fold cutoff. The default is
#' set to 2.
#' @param xlab \emph{x}-axis label.
#' @param ylab \emph{y}-axis label.
#' @param labelunderlim A logical indicating whether to label points that are not significant.
#' @param labelsig A logical indicating whether all significant points should be labeled.
#' @param interactiveplot A logical indication whether an interactive plot should be shown.
#' @param saveinteractiveplot A logical indication whether the interactive plot produced should
#' be saved as a \code{.html} file.
#' @param interactiveplotname A character string indicating the name to be used for saving the
#' interactive plot.
#' @param interactiveonly A boolean whether only an interactive version of the plot is required
#' @param main Plot title.
#' @param fclabel An optional character sting to label the vertical coefficient cutoff.
#' @param tolabel A list of metabolite names on the graph to be labeled
#' @param vlines A logical indicating whether to show vertical coefficient cutoff lines.
#' @param ... Other graphical parameters. See \code{\link[graphics]{par}}.
#' @param chooselegend Defualt to \code{NULL}. For internal use in other functions in the package.
#' @author Alysha M De Livera, Gavriel Olshansky
#' @examples
#'  data("alldata_eg")
#' logdata<-LogTransform(alldata_eg$featuredata)
#' sampledata<-alldata_eg$sampledata
#' metabolitedata<-alldata_eg$metabolitedata
#' imp <-  MissingValues(logdata$featuredata,sampledata,metabolitedata,
#'                      feature.cutof=0.8, sample.cutoff=0.8, method="knn")
#' featuredata<-imp$featuredata
#' qcmets<-which(metabolitedata[,1]=="IS")
#' factormat<-model.matrix(~gender +Age , sampledata)
#'
#' #Linear model fit with ordinary statistics with ruv2
#' ordFit_ruv2<-LinearModelFit(featuredata=featuredata,
#'                            factormat=factormat,
#'                            ruv2=TRUE, qcmets=qcmets,
#'                            k=2)
#' #Volcano plot
#' VolcanoPlot(coef=ordFit_ruv2$coefficients[,3], 
#'            pval=ordFit_ruv2$p.value[,3],
#'            cexlab = 0.8, 
#'            interactiveplot = TRUE, 
#'            coeflimit = 0.05,
#'            xlab="Coef",
#'            negcontrol= rownames(ordFit_ruv2$coefficients)
#'            [which(metabolitedata[,2]==1)],
#'            poscontrol= c("m74", "m161"),
#'            interactiveonly = TRUE)
#' @export VolcanoPlot
VolcanoPlot <- function(coef, pvals,
                        cexcutoff=0.7, 
                        cexlab=0.5,pointsize = 0.9, 
                        plimit=0.05, coeflimit=1, 
                        negcontrol = NULL, 
                        poscontrol=NULL,saveplot=FALSE, 
                        plotname="VolcanoPlot",
                        savetype= c("png","bmp","jpeg","tiff","pdf"),
                        xlab='Coefficients', ylab='-log(p-value)', 
                        labelunderlim = FALSE, labelsig=FALSE,
                        interactiveplot= TRUE, 
                        saveinteractiveplot=FALSE, 
                        interactiveplotname="interactiveVolcanPlot",
                        interactiveonly = FALSE,
                        main="Volcano Plot",fclabel="",
                        chooselegend = NULL,             #only pass if part of making comparison plot
                        vlines=TRUE, tolabel=NULL, ...){
  
  # Get range for x-vals
  #x_min <- (-1.5)
  #x_max <- 1.5
  #if (min(range(coef, finite = TRUE)) <= x_min) {
  #    x_min <- min(range(coef, finite = TRUE))
  #}
  #if (max(range(coef, finite = TRUE)) >= x_max) {
  #    x_max <- max(range(coef, finite = TRUE))
  #}
  
  x_min <- 1.2*min(c(coef,-1*coeflimit))
  x_max <- 1.2*max(c(coef,coeflimit))
  x_range <- c(x_min, x_max)
  
  
  # Get range for y-vals
  y_min <- 0
  y_max <- 2
  if (min(range(-log10(pvals), finite = TRUE)) <= y_min) {
    y_min <- min(range(-log10(pvals), finite = TRUE))
  }
  if (max(range(-log10(pvals), finite = TRUE)) >= y_max) {
    y_max <- max(range(-log10(pvals), finite = TRUE))
  }
  y_range <- c(y_min, y_max)
  
  if (class(coef) %in% c("data.frame", "list", "matrix")) 
    stop("coef should be a vector")

  if (class(pvals) %in% c("data.frame", "list", "matrix")) 
    stop("pvals should be a vector")
  
  #Stop if coefficient names are missing
  if(is.null(names(coef)))
    stop("coef vector must have names")
  
  
  #logical
  if(interactiveonly == TRUE){
    interactiveplot <- TRUE
  }
  ############################################
  ### set right legend label indicators and if variable passed for comparison plot
  ### Also sets parameter to know if part of multiplot
#  Defualt is \code{NULL}, otherwise pass a vector with 6 entries of \code{TRUE} or 
  #' \code{FALSE}, Indicating whether a legend should be drown for that item.
  if (is.null(chooselegend)){
    chooselegend <- rep(TRUE,6)
    formultiplot <- FALSE
  } else {
    formultiplot <- TRUE
  }
  ############################################
  
  
  # make a list to store the references for the different types of metabolites for plotly plot
  FCNS <- vector('numeric')
  NS <- vector('numeric')
  SIG <- vector('numeric')
  SIGFOLD <- vector('numeric')
  POSCONT <- vector('numeric')
  NEGCONT <- vector('numeric')
  
  
  if(interactiveonly == FALSE){
    
    # Start saving plot if needed
    if (saveplot == TRUE){
      savetype <- match.arg(savetype)
      savef <- match.fun(savetype)
      if (savetype == "pdf") {
        savef(paste(c(plotname,".",savetype),collapse = ""),width = 10, height = 9)
      }else {
        savef(paste(c(plotname,".",savetype),collapse = ""),width = 840, height = 720)
      }
    }
    
    # Draw the basic (empty) plot
    plot(
      x_range,                           # x-dim
      y_range,                           # y-dim
      type='n',                          # empty plot
      xlab=xlab,                         # x-axis title
      ylab=ylab,                         # y-axis title
      main=main,                         # plot title
      ...
    )
    
    # Annotate plot region
    abline(h=-log10(plimit),               # horizontal line at P=plimit
           col='green',                       # line colour
           lty='44'                           # Dot-dash lengths
    )
    mtext(paste("pval =", plimit),         # Label abline
          side=2,                            # on the left plot edge
          at=-log10(plimit),                 # at P=plimit
          cex=cexcutoff,                     # slightly smaller
          las=1                              # perpendicular to axis
    )
    # plot vertical line if user chooses to
    if (vlines == TRUE){
      abline(                                # vertical lines at plus minus 2-fold
        v=c(-(coeflimit), (coeflimit)),
        col='violet',
        lty='1343'
      )
      mtext(                                 # Label vertical ablines  #need to add option of label here
        c(paste("-", coeflimit, fclabel), paste("+", coeflimit, fclabel)),
        side=3,                            # on top of graph
        at=c(-coeflimit, coeflimit),
        cex=cexcutoff,
        las=1
      )
    }
    
    # write function that plots according to colour, and point (can save lots of space)
    # to replace the point() text() combination...
    
    # indicator whether to name the current point
    namepoint <- FALSE
    
    # Plot coloured points based on their values
    for (ii in 1:length(pvals)) {
      # Check if p-value is zero and if so add a small numebr to it
      
      # If part of our control groups, label differently
      if (is.element(names(coef)[ii],negcontrol)) {
        # update ref for interactive plot
        NEGCONT <- c(NEGCONT, ii)
        
        # Negative control so plot and colour:   green
        points(coef[ii],
               -log10(pvals[ii]),
               col='green',
               pch=20,cex= pointsize
        )
        
      } else if (is.element(names(coef)[ii],poscontrol)) {
        # update ref for interactive plot
        POSCONT <- c(POSCONT, ii)
        
        # Positive control so plot and colour:    red
        points(coef[ii],
               -log10(pvals[ii]),
               col='red',
               pch=20, cex= pointsize
        )
        
      } else {
        # If it's below plimit, we're not overly interested: purple.
        if (-log10(pvals[ii]) > (-log10(plimit))) {
          # Otherwise, more checks;
          # if it's greater than 2-fold decrease: blue
          if (coef[ii] > (-coeflimit)) {
            
            # If it's significant but didn't change much: orange
            if (coef[ii] < coeflimit) {
              # update ref for interactive plot
              SIG <- c(SIG, ii)
              
              points(coef[ii],
                     -log10(pvals[ii]),
                     col='orange',
                     pch=20,cex= pointsize
              )
              # points not very important, only label if user chooses so
              if (labelunderlim == TRUE){
                namepoint <- TRUE
              }
              # Otherwise, greater than 2-fold increase: blue
            } else {
              # update ref for interactive plot
              SIGFOLD <- c(SIGFOLD, ii)
              
              points(coef[ii],
                     -log10(pvals[ii]),
                     col='blue',
                     pch=20,cex = pointsize+0.3
              )
              if(labelsig == TRUE){
                namepoint<-TRUE
                
                
              }
            }
            # Else it's less than -2-fold decrease: blue
          } else {
            # update ref for interactive plot
            SIGFOLD <- c(SIGFOLD, ii)
            
            points(coef[ii],
                   -log10(pvals[ii]),
                   col='blue',
                   pch=20, cex = pointsize
            )
            if(labelsig == TRUE){
              namepoint <- TRUE
            }
          }
          # Else P > plimit; not significant: purple
        } else {
          # above coeflimit
          if(abs(coef[ii]) >= coeflimit){
            FCNS <- c(FCNS,ii)
            
            points(coef[ii],
                   -log10(pvals[ii]),
                   col='brown',
                   pch=20,cex = pointsize)
            # under coeflimit
          }else{
            # update ref for interactive plot
            NS <- c(NS,ii)
            
            points(coef[ii],
                   -log10(pvals[ii]),
                   col='purple',
                   pch=20,cex = pointsize
            )
          }
        }
      }
      if (is.element(names(coef)[ii],tolabel)) {
        namepoint <- TRUE
      }
      
      # only label points on graph that where choosen to
      if (namepoint==TRUE){
        text(coef[ii],        # x-coord
             -log10(pvals[ii]), # y-coord
             labels=names(coef)[ii],
             # If the point is at the top of the
             # graph, label goes underneath. If it's
             # at the far right, put the label on
             # the left of the point.
             
             pos=if(-log10(pvals[ii]) < 0.95 * max(y_range)) {
               if(coef[ii] < 0.75 * max(x_range)) {
                 4          # right if it's neither
               } else {
                 2          # left if > 0.75 max(x_range)
               }
             } else {
               1              # bottom if > 0.95 max(y_range)
             },
             cex=cexlab         # Size of text
        )
      }
      namepoint <-FALSE  
    }
    
    # Stop saving output
    if (saveplot==TRUE){
      dev.off()
    }
    
    # only need an interactive plot so get the right grouping for the metabolites 
  } else{                                  
    
    for (ii in 1:length(pvals)) {
      # Check if p-value is zero and if so add a small numebr to it
      
      # If part of our control groups, label differently
      if (is.element(names(coef)[ii],negcontrol)) {
        # update ref for interactive plot
        NEGCONT <- c(NEGCONT, ii)
        
      } else if (is.element(names(coef)[ii],poscontrol)) {
        # update ref for interactive plot
        POSCONT <- c(POSCONT, ii)
        
        
      } else {
        # If it's below plimit, we're not overly interested: purple.
        if (-log10(pvals[ii]) > (-log10(plimit))) {
          # Otherwise, more checks;
          # if it's greater than 2-fold decrease: blue
          if (coef[ii] > (-coeflimit)) {
            
            # If it's significant but didn't change much: orange
            if (coef[ii] < coeflimit) {
              # update ref for interactive plot
              SIG <- c(SIG, ii)
              
              
              # Otherwise, greater than 2-fold increase: blue
            } else {
              # update ref for interactive plot
              SIGFOLD <- c(SIGFOLD, ii)
              
            }
            # Else it's less than -2-fold decrease: blue
          } else {
            # update ref for interactive plot
            SIGFOLD <- c(SIGFOLD, ii)
            
          }
          # Else P > plimit; not significant: purple
        } else {
          # above coeflimit
          if(abs(coef[ii]) >= coeflimit){
            FCNS <- c(FCNS,ii)
            
            # under coeflimit
          }else{
            # update ref for interactive plot
            NS <- c(NS,ii)
            
          }
        }
      }
    }
  }
  
  if (interactiveplot == TRUE){
    if (interactiveonly == FALSE){
      cat("Hit <Return> to see interactive plot: ")
      line <- readline()
    }
    
    text1 <- if(length(POSCONT)==0){
      NULL
    } else {
      paste("Metabolite: ",names(coef)[POSCONT],
            "<br> P-Value: ",round(pvals[POSCONT],3),     #set sig figs for p-Value
            "<br> Coefficient: ",round(coef[POSCONT],2),  #set sig figs for coef
            sep = "")
    }
    text2 <- if(length(NEGCONT)==0){
      NULL
    } else { 
      paste("Metabolite: ",names(coef)[NEGCONT],
            "<br> P-Value: ",round(pvals[NEGCONT],3),
            "<br> Coefficient: ",round(coef[NEGCONT],2),
            sep = "")
    }
    text3 <- if(length(NS)==0){
      NULL
    } else {
      paste("Metabolite: ",names(coef)[NS],
            "<br> P-Value: ",round(pvals[NS],3),
            "<br> Coefficient: ",round(coef[NS],2),
            sep = "")
    }
    text4 <- if(length(SIG)==0){
      NULL
    } else {
      paste('Metabolite: ',names(coef)[SIG],
            '<br> P-Value: ',round(pvals[SIG],3),
            '<br> Coefficient: ',round(coef[SIG],2)
            )
    }
    text5 <- if(length(SIGFOLD)==0){
      NULL
    } else {
      paste("Metabolite: ",names(coef)[SIGFOLD],
            "<br> P-Value: ",round(pvals[SIGFOLD],3),
            "<br> Coefficient: ",round(coef[SIGFOLD],2),
            sep = "")
    }
    text6 <- if(length(FCNS)==0){
      NULL
    } else {
      paste("Metabolite: ",names(coef)[FCNS],
            "<br> P-Value: ",round(pvals[FCNS],3),
            "<br> Coefficient: ",round(coef[FCNS],2),
            sep = "")
    }
    
    # positive control
    trace1 <- list(
      x = coef[POSCONT],
      y = c(-log10(pvals[POSCONT])),
      marker = list(color = "red"), 
      mode = "markers",
      name = "Positive Control",
      hoverinfo = "text",
      text = text1,
      textposition = "bottom",
      type = "scatter"
      #  uid = "3f9a46"
    )
    # negative control
    trace2 <- list(
      x = coef[NEGCONT],
      y = c(-log10(pvals[NEGCONT])),
      marker = list(color = "green"), 
      mode = "markers",
      name = "Negative Control",
      hoverinfo = "text",
      text = text2,
      textposition = "bottom",
      type = "scatter"
      #  uid = "c2d206"
    )
    # not significant results
    trace3 <- list(
      x = coef[NS],
      y = c(-log10(pvals[NS])),
      marker = list(color = "purple"), 
      mode = "markers",
      name = paste("p-value > ", plimit, "&", "|", xlab, "|", "<", coeflimit), #changed ADL 
      hoverinfo = "text",
      text = text3,
      textposition = "bottom",
      type = "scatter"
      #  uid = "37b094"
    )
    # significant result but low coefficient
    trace4 <- list(
      x = coef[SIG],
      y = c(-log10(pvals[SIG])),
      marker = list(color = "orange"), 
      mode = "markers",
      name = paste("p-value < ", plimit, "&", "|", xlab, "|", "<", coeflimit),  #Changed ADL 
      hoverinfo = "text",
      text = text4,
      textposition = "bottom",
      type = "scatter"
      #  uid = "37b094"
    )
    # significant results and coefficients
    if(labelsig == TRUE){
      sigmode = "markers+text"
    } else{
      sigmode= "markers"
    }
    
    trace5 <- list(
      x = coef[SIGFOLD],
      y = c(-log10(pvals[SIGFOLD])),
      marker = list(color = "blue"),
      mode = sigmode,
      name = paste("p-value < ", plimit, "&", "|", xlab, "|", ">", coeflimit), #changed ADL  
      hoverinfo = "text",
      text = text5,
      textposition = "bottom",
      type = "scatter"
      #  uid = "37b094"
    )
    # coefficients above min but not significant
    trace6 <- list(
      x = coef[FCNS],
      y = c(-log10(pvals[FCNS])),
      marker = list(color = "brown"), 
      mode = "markers",
      name = paste("p-value > ", plimit, "&", "|", xlab, "|", ">", coeflimit), #changed ADL 
      hoverinfo = "text",
      text = text6,
      textposition = "bottom",
      type = "scatter"
      #  uid = "37b094"
    )
    # add vertical lines if required
    if (vlines==TRUE){
      line1 <- list(
        type=line
      )
    }
    
    layout <- list(
      autosize = TRUE,
      showlegend = TRUE,
      title = main,
      #height = 602,   Set size to correspond to page size
      #width = 936,
      xaxis = list(
        anchor = "y",
        autorange = TRUE,
        autostick = TRUE,
        side = "bottom",
        title = xlab, #Changed adl 03/03/2017
        type = "linear"
      ),
      yaxis = list(
        autorange = TRUE,
        title = ylab, #"log(P-Value)", Changed ADL
        type = "linear"
      )
    )
    p<-plot_ly()
   
    if(length(NS)>0){            # add trace if have non significant results
      p <- add_trace(p, x=trace3$x, y=trace3$y, marker=trace3$marker, mode=trace3$mode, name=trace3$name, 
                   hoverinfo=trace3$hoverinfo, text=trace3$text, textposition=trace3$textposition, type=trace3$type,
                   legendgroup = "groupNS", showlegend = chooselegend[3])
      chooselegend[3] <- FALSE      # update incase making another plot in compare function
    }
    if(length(SIG)>0){       # add trace if have significant results with low coefficients
      p <- add_trace(p, x=trace4$x, y=trace4$y, marker=trace4$marker, mode=trace4$mode, name=trace4$name,
                   hoverinfo=trace4$hoverinfo, text=~trace4$text, textposition=trace4$textposition, type=trace4$type,
                   legendgroup = "groupSIG", showlegend = chooselegend[4])
      chooselegend[4] <- FALSE      # update incase making another plot in compare function
    }
    if(length(SIGFOLD)>0){   # add trace if have significant results with coefficients above..
      p <- add_trace(p, x=trace5$x, y=trace5$y, marker=trace5$marker, mode=trace5$mode, name=trace5$name,
                   hoverinfo=trace5$hoverinfo, text=trace5$text, textposition=trace5$textposition, type=trace5$type,
                   legendgroup = "groupSIGFOLD", showlegend = chooselegend[5])
      chooselegend[5] <- FALSE
    }
    if(length(FCNS)>0){   # add trace if have non significant results with high coeficients
      p <- add_trace(p, x=trace6$x, y=trace6$y, marker=trace6$marker, mode=trace6$mode, name=trace6$name,
                   hoverinfo=trace6$hoverinfo, text=trace6$text, textposition=trace6$textposition, type=trace6$type,
                   legendgroup = "groupFCNS", showlegend = chooselegend[6])
      chooselegend[6] <- FALSE
    }
    if(length(POSCONT)>0){       # add trace only if have positive controls
      p <- add_trace(p, x=trace1$x, y=trace1$y, marker=trace1$marker, mode=trace1$mode, name=trace1$name,
                     hoverinfo=trace1$hoverinfo, text=trace1$text, textposition=trace1$textposition, type=trace1$type,
                     legendgroup = "groupPOSCONT", showlegend = chooselegend[1])
      chooselegend[1] <- FALSE
    }
    if(length(NEGCONT)>0){       # add trace if have negative controls
      p <- add_trace(p, x=trace2$x, y=trace2$y, marker=trace2$marker, mode=trace2$mode, name=trace2$name,
                     hoverinfo=trace2$hoverinfo, text=trace2$text, textposition=trace2$textposition, type=trace2$type,
                     legendgroup = "groupNEGCONT", showlegend = chooselegend[2])
      chooselegend[2] <- FALSE
    }
    
    #p <- add_trace(p, x = c(coeflimit,coeflimit),x =c(-coeflimit,-coeflimit), y = c(0,-log10(min(pvals))),y = c(0,-log10(min(pvals))),name="Coefficient cutoff", type="scatter",mode="lines",line = list(dash = "dash"))
    p <- layout(p, autosize=layout$autosize, showlegend=layout$showlegend, 
                title=layout$title, xaxis=layout$xaxis, yaxis=layout$yaxis)
    
    if (saveinteractiveplot){                                               #save interactive plot if needed
      htmlwidgets::saveWidget(p, paste(interactiveplotname,".html",sep="")) #Changed ADL 03/03/2017
    }
    if (formultiplot==FALSE){
      return(p)
    } else {
      return(list(p,chooselegend))   ##need the extra parameter to pass to next plot
    }
  }
  
}













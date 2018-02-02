#' PCA plots
#' 
#' Produces PCA plots of the metabolomics data.
#' 
#' 
#' @param featuredata featuredata A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' See NormalizeMets Vignette for details.
#' @param groupdata A data frame or a vector with group names.
#' @param saveplot A logical indication whether to save the plot produced.
#' @param saveinteractiveplot A logical indication whether to save the interactive plots produced.
#' @param plotname Name of the output file if the file is to be saved. This is
#' the general name for all the graphs and the specific type prefix will be
#' added automatically.
#' @param savetype The required format for the plot to be saved in. Threre is a
#' choice of \code{"png","bmp","jpeg","tiff","pdf"} type files.
#' @param interactiveonly Alogical indicating whether to show interactive 
#' plots only.
#' @param interactiveplots A logical indication whether an interactive plot
#' is to be shown.
#' @param y.axis The principal component to be plotted on the \emph{y}-axis.
#' @param x.axis The principal component to be plotted on the \emph{x}-axis.
#' @param center A logical indicating whether the variables should be scaled to
#' have zero mean.
#' @param scale A logical indicating whether the variables should be scaled to
#' have unit variance before the analysis takes place.
#' @param userinput A logical indication whether user input should be required to 
#' show plots progressively. Should be \code{FALSE} when the plots are to be passed
#' on to other functions or saved as variables for later use.
#' @param returninteractive A logical indication whether a list of the interactive plots
#' should be returned by the function
#' @param main Plot title.
#' @param varplot A logical indicating whether explained variance should be
#' plotted.
#' @param multiplot If \code{TRUE}, pairs plots of the first \emph{n} principal
#' components will be plotted.
#' @param n The number of principal components to be plotted if
#' \code{multiplot=TRUE}. The default value is set to 5.
#' @param cols A character string with colours to be used.
#' @param cex_val A numeric indicating the size of some text elements.
#' @param ... Arguments to be passed on to other methods.
#' @author Alysha M De Livera, Gavriel Olshansky
#' @seealso \code{\link[stats]{prcomp}}.
#' @examples
#' 
#'     data(mixdata)
#'     
#'     # produce all results
#'     PcaPlots(mixdata$featuredata,mixdata$sampledata[,3],multiplot = TRUE, 
#'     varplot = TRUE, interactiveplots = TRUE)
#'     
#'     # return a list of the ineractive plots only
#'     interactive.pca <- PcaPlots(mixdata$featuredata,mixdata$sampledata[,3],
#'     interativeonly = TRUE, interactiveplots = TRUE,
#'     userinput = FALSE, returninteractive = TRUE)
#'     
#' 
#' @export PcaPlots
PcaPlots <- function(featuredata, groupdata, saveplot=FALSE,saveinteractiveplot = FALSE, 
                     plotname="",savetype= c("png","bmp","jpeg","tiff","pdf"), interactiveonly = FALSE,
                     interactiveplots = TRUE, y.axis=1, x.axis=2, center=TRUE, scale=TRUE,
                     userinput = TRUE, returninteractive = FALSE,
                     main=NULL, varplot=FALSE, multiplot=FALSE, n=3, cols=NULL,cex_val = 0.7, ...)
{
  
  # prepare plot names for saving
  if (saveplot == TRUE | saveinteractiveplot == TRUE){
    savetype <- match.arg(savetype)
    plottype <- c("varplot","multiplot","PCA_score","PCA_loading")
    savenames <- vector( ,4)
    
    # edit name for saving
    if(length(plotname) != 0) {
      plotname <- paste(plotname,"_",sep = "")
    }
    for (i in 1:4){
      savenames[i] <- paste(plotname,plottype[i],".",savetype,sep = "")
    }
  }
  
  # Get groups information
  group_list <- factor(groupdata, levels=unique(groupdata))
  # Remove groups for data processing
  pca_data <- featuredata
  
  if(userinput){
    write(' -> Performing PCA...', '')
  }
  
  const_rows <- which(apply(pca_data, 2, var) == 0)
  if (length(const_rows) != 0) {
    pca_data <- pca_data[, -const_rows]
  }
  
  pca <- prcomp(pca_data, scale.=scale, center=center, ...)
  # Get the eigenvectors (loadings)
  eigenvecs <- pca$rotation
  
  # Get summary information
  summ <- summary(pca)
  importance <- summ$importance[2, ]
  
  
  if (varplot) {
    # Plot the explained variance
    # dev.new()
    # Save if required
    if (saveplot == TRUE) {
      savef <- match.fun(savetype)
      if (savetype == "pdf") {
        savef(savenames[1],width = 10, height = 9)
      }else {
        savef(savenames[1],width = 840, height = 720)
      }
    }
    # MAke the var Plot
    barplot(summ$importance[2, c(1:10)],
            col="#ee3333",                  # colour to plot bars
            main="Variance",               # plot title
            xlab="Principal Component",    # x-axis title
            ylab="Explained Variance",     # y-axis title
            cex.axis= 1,
            cex.names=1,       # font size
            las=1                          # horizontal labels
    )
    # Stop saving to current file
    if (saveplot == TRUE) { 
      dev.off()
    }
  }
  
  # Prepare a list of colours to use
  if (is.null(cols)) {
    col_list <- ColList(nlevels(group_list))
    cols_used <- col_list[as.numeric(group_list)]
  } else
    cols_used<-cols[as.numeric(group_list)]
  pch_list <- PchList(nlevels(group_list))
  pch_used <- pch_list[as.numeric(group_list)]
  
  # Get info for ploting legend
  unique.groups <- levels(group_list)
  cols <- ColList(length(unique.groups))
  
  # Store PCA data in a meaningful namespace
  pca_scores <- pca$x
  rownames(pca_scores) <- rownames(pca_data)
  
  
  
  #Create the Multiplot 
  
  if (multiplot) {
    
    #check if plot needs to be saved
    if (saveplot == TRUE) {
      if (savetype == "pdf") {
        savef(savenames[2],width = 10, height = 9)
      }else {
        savef(savenames[2],width = 840, height = 720)
      }
      #wait for key stroke to show next plot
    } else if (userinput) {
      cat("Hit <Return> to see next plot: ")
      line <- readline()
    }
    
    #dev.new()
    pairs(pca_scores[,1:n], pch=pch_used, col=cols_used,
          labels=paste(
            "PC", c(1:n), "(", round(importance[c(1:n)] * 100, 2), "%)",
            sep="")
    )
    
    # Stop saving to current file
    if (saveplot == TRUE) { 
      dev.off()
    }
  }
  
  # Plot PCA scores with sample names
  pca_mat <- cbind(pca_scores[, x.axis],pca_scores[, y.axis])
  x_percent <- sprintf("%.2f", importance[x.axis] * 100)
  y_percent <- sprintf("%.2f", importance[y.axis] * 100)
  
  
  #check if plot needs to be saved
  if(interactiveonly == FALSE){
    if (saveplot == TRUE) {
      if (savetype == "pdf") {
        savef(savenames[3],width = 10, height = 9)
      }else {
        savef(savenames[3],width = 840, height = 720)
      }
      #wait for key stroke to show next plot
    } else if (userinput) {
      cat("Hit <Return> to see next plot: ")
      line <- readline()
    }
    
    pic_gen(pca_mat,
            plot_title=if (!is.null(main)) main else "PCA Score Plot\nSamples",
            x_label= paste("PC", x.axis, " (", x_percent, "%)", sep=""),
            y_label=paste("PC", y.axis, " (", y_percent, "%)", sep=""),
            cols_used=cols_used,
            pch_used=pch_used,cex_val=cex_val
    )
    # Add legend to the graph
    legend("bottomleft", legend = unique.groups,
           col=cols, lty = c(0,0),pch=c(15,19),lwd = c(1.5,1.5),cex = cex_val)
    
    # Stop saving to current file
    if (saveplot == TRUE) { 
      dev.off()
    }
  }
  
  
  # Not needed as a plot with sample names has the group data as well (by colour)
  # Plot PCA scores with group names
  #pic_gen(pca_mat,
  #   plot_title=if (!is.null(main)) main else "PCA Score Plot\nGroups",
  #    plot_labels=group_list,
  #    x_label=paste("PC", x.axis, " (", x_percent, "%)", sep=""),
  #    y_label=paste("PC", y.axis," (", y_percent, "%)",  sep=""),
  #    cols_used=cols_used,
  #    pch_used=pch_used
  #)
  
  
  eigen_mat <- cbind(eigenvecs[, x.axis],eigenvecs[, y.axis])
  
  if (interactiveonly== FALSE){
    
    #check if plot needs to be saved
    
    if (saveplot == TRUE) {
      if (savetype == "pdf") {
        savef(savenames[4],width = 10, height = 9)
      }else {
        savef(savenames[4],width = 840, height = 720)
      }
      #wait for key stroke to show next plot
    } else if (userinput){
      cat("Hit <Return> to see next plot: ")
      line <- readline()
    }
    
    
    # Loadings plot
    pic_gen(eigen_mat,
            "PCA Loading Plot",
            x_label=paste("PC", x.axis," (", x_percent, "%)", sep=""),
            y_label=paste("PC", y.axis," (", y_percent, "%)", sep=""),
            cols_used="black", cex_val=cex_val
            #text_on=FALSE
    )
    # Stop saving to current file
    if (saveplot == TRUE) { 
      dev.off()
    }
  }
  
  # Make interactive plots if required
  if(interactiveplots == TRUE){
    
    hover_infoS <- c()
    hover_infoL <- c()
    
    for (ii in 1:length(pca_mat[,1])){
      #score
      hover_infoS[ii] <- paste0("PC", x.axis,": ",pca_mat[ii,1],
                              "<br>PC", y.axis,": ",pca_mat[ii,2],
                              "<br>Sample number: ", rownames(pca_mat)[ii])
    }
    for (ii in 1:length(eigen_mat[,1])){
      #loading
      hover_infoL[ii] <- paste0("Loading PC", x.axis,": ",eigen_mat[ii,1],
                               "<br>Loading PC", y.axis,": ",eigen_mat[ii,2],
                               "<br>Metabolite number: ", rownames(eigen_mat)[ii])
    }
    
    gtoname<- setNames(unique(cols_used),c(unique.groups))
    
    i_score <- plot_ly(x = ~pca_mat[,1],y= ~pca_mat[,2],color = ~group_list, colors = gtoname,
                       hoverinfo = "text",
                       text = ~hover_infoS,
                       type = "scatter", mode = "markers") %>%
      layout(title = "PCA Score Plot",
             xaxis = list(title = paste("PC", x.axis, " (", x_percent, "%)", sep="")),
             yaxis = list(title = paste("PC", y.axis, " (", y_percent, "%)", sep="")))
    #plot score plot
    if (userinput){
      cat("Hit <Return> to see interactive plot: ")
      line <- readline()
    }
    # only show plot if not returning it
    if(!returninteractive){
      print(i_score)
    }
    
    i_loading <- plot_ly(x = ~eigen_mat[,1],y= ~eigen_mat[,2],
                         hoverinfo = "text",
                         text = ~hover_infoL,
                         type = "scatter",mode = "markers") %>%
      layout(title = "PCA Loading Plot",
             xaxis = list(title = paste("PC", x.axis, " (", x_percent, "%)", sep="")),
             yaxis = list(title = paste("PC", y.axis, " (", y_percent, "%)", sep="")))
    
    # plot loading plot
    
    if (userinput){
      cat("Hit <Return> to see next interactive plot: ")
      line <- readline()
    }
    if (saveinteractiveplot){                                                   #GO 29/5/17
      htmlwidgets::saveWidget(i_score, paste(plotname,plottype[3],".html",sep=""),selfcontained = F)
      htmlwidgets::saveWidget(i_loading, paste(plotname,plottype[4],".html",sep=""),selfcontained = F)
    }                  #Need selfcontain = FALSE to save both plots when running in shell
    
    # only show plot if not returning it
    if(!returninteractive){
      print(i_loading)
    }
  }
  if (userinput){
    write(' -> Done!', '')
  }
  
  if(returninteractive == TRUE){
    return(list(i_score, i_loading))
  }
  
  
  
}
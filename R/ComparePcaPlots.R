#' Compare PCA Plots
#' 
#' Produces a comparison of principal component multiplots
#' 
#' 
#' @param lfeaturedata A list of data frames in the featuredata format. 
#' This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' See NormalizeMets Vignette for details.
#' @param lgroupdata A list of data frames or a vectors with group names.
#' @param saveplot A logical indication whether to save the plot produced.
#' @param plotname Name of the output file if the file is to be saved. This is
#' the general name for all the graphs and the specific type prefix will be
#' added automatically.
#' @param savetype The required format for the plot to be saved in. Threre is a
#' choice of \code{"png","bmp","jpeg","tiff","pdf"} type files.
#' @param y.axis The principal component to be plotted on the \emph{y}-axis.
#' @param x.axis The principal component to be plotted on the \emph{x}-axis.
#' @param center A logical indicating whether the variables should be scaled to
#' have zero mean.
#' @param scale A logical indicating whether the variables should be scaled to
#' have unit variance before the analysis takes place.
#' @param lmain A list of plot titles
#' @param n The number of principal components to be plotted. The default value 
#' is set to 3.
#' @param showlegend A logical indication whether to print a legend for the plot.
#' @param usercols A character string with colours to be used or TRUE for ColList
#' to be used (ColList is the defualt set of colours used in other plots in this 
#' package). A value other then \code{NULL} will automatically set \code{showlegend} to 
#' \code{FALSE}.
#' @param cex_val A numeric indicating the size of some text elements.
#' @param ... Arguments to be passed on to other methods.
#' @author Alysha M De Livera, Gavriel Olshansky.
#' 
#' 
#' @examples
#' 
#'     data(mixdata)
#' #     ComparePcaPlots(list(mixdata$featuredata,mixdata$featuredata*1.2),
#' #             list(mixdata$sampledata[,3],mixdata$sampledata[,3])) 
#' 
ComparePcaPlots <- function(lfeaturedata, lgroupdata, saveplot=FALSE, 
                     plotname="",savetype= c("png","bmp","jpeg","tiff","pdf"),
                     y.axis=1, x.axis=2, center=TRUE, scale=TRUE,
                     lmain=NULL, n=3, showlegend =TRUE, usercols=NULL,
                     cex_val = 0.7, ...)
{
  
  # prepare plot names for saving
  if (saveplot == TRUE ){
    savetype <- match.arg(savetype)
    plottype <- c("multiplot")     #"varplot","PCA_score","PCA_loading")
    savenames <- vector( ,1)
    
    #edit name for saving
    #if(length(plotname) != 0) {
    # plotname <- paste(plotname,"_",sep = "")
    #}
    for (i in 1:1){   # for now
      savenames[i] <- paste(plotname,"_",plottype[i],".",savetype,sep = "")
    }
  }
  
  # Get required info for multiple plots #
  
  num_of_plots <- length(lfeaturedata)
  
  # if only one grouping entered, use it for all plots
  if (length(lgroupdata) ==1 & (num_of_plots>1)){
    my.group <- lgroupdata[[1]]
    for (ii in 1:num_of_plots){
      lgroupdata[[ii]] <- my.group
    }
  }
  
  if (length(lgroupdata)!= num_of_plots){
    stop("Make sure that lgroupdata is of the same size as featuredatal")
  }
  if (num_of_plots > 4){
    stop("Too many plots! (4 is the maximum for comparison)")
  }
  if (length(lmain) < num_of_plots){
    lmain <- list()
    for (ii in 1:num_of_plots){
      lmain[[ii]] <- paste("Method_",ii,sep="")
    }
  }
  
  plotnumcol <- ifelse(num_of_plots > 1,2,1)  # number of columns to use for plots
  plotnumrow <- ifelse(num_of_plots > 2,2,1)
  
 
  # prepare variables to save info for multiple plots
  group_list <- list()
  pca_data <- list()
  pca <- list()
  eigenvecs <- list()
  summ <- list()
  importance <- list()
  
  write(' -> Performing PCA...', '')
  
  # Get info and factors/groups for plots - do the pc breakdown
  for (jj in 1:num_of_plots){
  
    # Get groups information
    group_list[[jj]] <- factor(lgroupdata[[jj]], levels=unique(lgroupdata[[jj]]))
    # Remove groups for data processing
    pca_data[[jj]] <- lfeaturedata[[jj]]
  
    const_rows <- which(apply(pca_data[[jj]], 2, var) == 0)
    if (length(const_rows) != 0) {
      pca_data[[jj]] <- pca_data[[jj]][, -const_rows]
    }
  
    pca[[jj]] <- prcomp(pca_data[[jj]], scale.=scale, center=center, ...)
    # Get the eigenvectors (loadings)
    eigenvecs[[jj]] <- pca[[jj]]$rotation
  
    # Get summary information
    summ[[jj]] <- summary(pca[[jj]])
    importance[[jj]] <- summ[[jj]]$importance[2, ]
  }
  
  
  #if (varplot) {
  # # Plot the explained variance
  #  # dev.new()
  #  # Save if required
    
  #  if (saveplot == TRUE) {
  #    savef <- match.fun(savetype)
  #    if (savetype == "pdf") {
  #      savef(savenames[1],width = 10, height = 9)
  #    }else {
  #      savef(savenames[1],width = 840, height = 720)
  #    }
  #  }
  #  ### Make the var Plot ###
    
    #save and change par settings for multiple plots
  #  oldpar <- par()
  #  par(mfrow = mymfrow)
    
  #  for (jj in 1:num_of_plots){
  #    barplot(summ[[jj]]$importance[2, c(1:10)],
  #            col="#ee3333",                  # colour to plot bars
  #            main= paste("Variance ",lmain[[jj]],sep=""),         # plot title
  #            xlab="Principal Component",    # x-axis title
  #            ylab="Explained Variance",     # y-axis title
  #            cex.axis= 1,
  #            cex.names=1,       # font size
  #            las=1                          # horizontal labels
  #    )
  #  }
  #  par(oldpar)         # reset par settings to par before plot
    
  #  # Stop saving to current file
  #  if (saveplot == TRUE) { 
  #    dev.off()
  #  }
  #}
  ################
  
  
  # Prepare a list of colours to use
  
  col_list <- list()
  cols_used <- list()
  pch_list <- list()
  pch_used <- list()
  unique.groups <- list()
  cols <- list()
  pca_scores <- list()
  pca_gg <- list()   # used for multiplot with ggpairs
  max.cols <- 0  # record the max colours needed
  second.pos <- c(5:7) # the defualt position for the second plot in a row
  my.leg.loc <- c(1,3)
  
  for (jj in 1:num_of_plots){
    
    col_list[[jj]] <- ColList(nlevels(group_list[[jj]]))  ##ColList
    cols_used[[jj]] <- col_list[[jj]][as.numeric(group_list[[jj]])]
  
    pch_list[[jj]] <- PchList(nlevels(group_list[[jj]]))
    pch_used[[jj]] <- pch_list[[jj]][as.numeric(group_list[[jj]])]

    # Get info for ploting legend
    unique.groups[[jj]] <- levels(group_list[[jj]])
    cols[[jj]] <- ColList(length(unique.groups[[jj]]))   #ColList
  
    # Store PCA data in a meaningful namespace
    pca_scores[[jj]] <- pca[[jj]]$x
    rownames(pca_scores[[jj]]) <- rownames(pca_data[[jj]])
    Grouping <- list()   # define an empty list to store group information in functions
    # for multiplot
    pca_gg[[jj]] <- cbind(as.data.frame(pca_scores[[jj]]),Grouping = lgroupdata[[jj]])
    
    # update max colours if needed
    max.cols <- max(max.cols, length(unique.groups[[jj]]))
  }
  
  ###### Check if need to use custom colours #######
  
  if (!is.null(usercols)){
    
    showlegend <- FALSE
    UserCols <- TRUE
    my.leg.loc <- NULL
    
    if (usercols == TRUE){
      CustomCols <- ColList(max.cols)
    } else {
      if (length(usercols) < max.cols){
        stop("Not enough custom colours given!")
      } else{
        CustomCols <- usercols
      }
    }
  }
  
  if (!showlegend)
    second.pos <- c(4:6)
  
  
  ###### Create the Multiplot ####### 
  multiplot <-  TRUE     # for now always make a multiplot
  
  if (multiplot) {
    
    pairgg <- list()            # use this for the individual plots
    
    #check if plot needs to be saved
    if (saveplot == TRUE) {
      savef <- match.fun(savetype)
      if (savetype == "pdf") {
        savef(savenames[1],width = 6*plotnumcol+2, height = 6*plotnumrow)
      }else {
        savef(savenames[1],width = 480*plotnumcol+160, height = 480*plotnumrow)
      }
      #wait for key stroke to show next plot
    }
    
    ##### Try adding a plot made with ggpairs ######
    if (num_of_plots == 1){
      #only one plot so can have legend directly
      pairgg[[1]] <- ggpairs(pca_gg[[jj]],columns = 1:n,aes(colour = Grouping),
                             upper = list(continuous="points"), legend = my.leg.loc,
                             title = lmain[[jj]]) + ggplot2::theme(legend.position = "right")
              
      if (!is.null(usercols)){      #add custom cols if needed
        for(i in 1:n){
          for(j in 1:n){
            pairgg[[1]][i,j] <- pairgg[[1]][i,j] + scale_color_manual(values = CustomCols) + scale_fill_manual(values = CustomCols)
          }
        }
      }
      
      # get top left y-scale and bottom rigth x scale for use in plot
      my.yrange <- c(min(ggplot_build(pairgg[[jj]][1,n])$data[[1]]$y),max(ggplot_build(pairgg[[jj]][1,n])$data[[1]]$y))
      my.xrange <- c(min(ggplot_build(pairgg[[jj]][1,n])$data[[1]]$x),max(ggplot_build(pairgg[[jj]][1,n])$data[[1]]$x))
      
      
      # add % var explained
      for (ii in 1:n){
        textgg <- ggally_text(paste(
          "PC", ii , "(", round(importance[[jj]][ii] * 100, 2), "%)",
          sep=""), xrange = my.xrange, yrange = my.yrange)
        pairgg[[1]][ii,ii] <- textgg
      }
      
      return(pairgg[[1]])
      
    }   #if 1 plot, return the one plot with legend
    
    ######################
    # if more than one plot, set the legend in the middle and all the other plots around it
    # first just generate plots without legend
    for (jj in 1:num_of_plots){
      
      pairgg[[jj]] <- ggpairs(pca_gg[[jj]],columns = 1:n,aes(colour = Grouping),
                              upper = list(continuous="points"),
                              title = lmain[[jj]])
                              
      if (!is.null(usercols)){      #add custom cols if needed
        for(i in 1:n){
          for(j in 1:n){
            pairgg[[jj]][i,j] <- pairgg[[jj]][i,j] + scale_color_manual(values = CustomCols) + scale_fill_manual(values = CustomCols)
          }
        }
      }
      
      # get top left y-scale and bottom rigth x scale for use in plot
      my.yrange <- c(min(ggplot_build(pairgg[[jj]][1,n])$data[[1]]$y),max(ggplot_build(pairgg[[jj]][1,n])$data[[1]]$y))
      my.xrange <- c(min(ggplot_build(pairgg[[jj]][1,n])$data[[1]]$x),max(ggplot_build(pairgg[[jj]][1,n])$data[[1]]$x))
      
      # add % var explained
      for (ii in 1:n){
        textgg <- ggally_text(paste(
          "PC", ii , "(", round(importance[[jj]][ii] * 100, 2), "%)",
          sep=""), xrange = my.xrange, yrange = my.yrange)
        pairgg[[jj]][ii,ii] <- textgg
      }
    }
    
    #### make legend for all plots #####
    points_legend <- gglegend(ggally_points)
    my_legend <- points_legend(pca_gg[[1]],mapping =ggplot2::aes("PC1","PC2",colour = Grouping))
    
    ###make plots with the correct layout #####
    if (num_of_plots == 2){
      grid.newpage()
      if (showlegend == TRUE){
        pushViewport(viewport(layout = grid.layout(1,7)))
      } else {
        pushViewport(viewport(layout = grid.layout(1,6)))
      }
      if (showlegend) 
        print(my_legend)
      print(pairgg[[1]], vp = viewport(layout.pos.row = 1, layout.pos.col = c(1:3)))
      print(pairgg[[2]], vp = viewport(layout.pos.row = 1, layout.pos.col = second.pos))
    } else {
      grid.newpage()
      if (showlegend == TRUE){
        pushViewport(viewport(layout = grid.layout(2,7)))
      } else {
        pushViewport(viewport(layout = grid.layout(2,6)))
      }
      if (showlegend) 
        print(my_legend)
      print(pairgg[[1]], vp = viewport(layout.pos.row = 1, layout.pos.col = c(1:3)))
      print(pairgg[[2]], vp = viewport(layout.pos.row = 1, layout.pos.col = second.pos))
      print(pairgg[[3]], vp = viewport(layout.pos.row = 2, layout.pos.col = c(1:3)))
      if (num_of_plots == 4){ 
        print(pairgg[[4]], vp = viewport(layout.pos.row = 2, layout.pos.col = second.pos))
      }
    }
    
    # Stop saving to current file
    if (saveplot == TRUE) { 
      dev.off()
    }
  }
  #############################################
  
  ## Plot PCA scores with sample names
  #pca_mat <- list() 
  #x_percent <- list()
  #y_percent <- list()
  
  # set parameters for next plot
  #for (jj in 1:num_of_plots){
  #pca_mat[[jj]] <- cbind(pca_scores[[jj]][, x.axis],pca_scores[[jj]][, y.axis])
  #x_percent[[jj]] <- sprintf("%.2f", importance[[jj]][x.axis] * 100)
  #y_percent[[jj]] <- sprintf("%.2f", importance[[jj]][y.axis] * 100)
  #}
  
  ######### Make the pca scores plot #############
  
  #check if plot needs to be saved
  #if (saveplot == TRUE) {
  #  if (savetype == "pdf") {
  #    savef(savenames[3],width = 10, height = 9)
  #  }else {
  #    savef(savenames[3],width = 840, height = 720)
  #  }
  #  #wait for key stroke to show next plot
  #} else {
  #  cat("Hit <Return> to see next plot: ")
  #  line <- readline()
  #}
  
  #oldpar <- par()           #set par settings for the plot
  #par(mfrow = mymfrow)
  
  ## make the plots
  #for (jj in 1:num_of_plots){
  #  pic_gen(pca_mat[[jj]],
  #          plot_title= paste("PCA Score Plot\nSamples - ",lmain[[jj]],sep=""),
  #          x_label= paste("PC", x.axis, " (", x_percent[[jj]], "%)", sep=""),
  #          y_label=paste("PC", y.axis, " (", y_percent[[jj]], "%)", sep=""),
  #          cols_used=cols_used[[jj]],
  #          pch_used=pch_used[[jj]],cex_val=cex_val
  #  )
  #  # Add legend to the graph -- Need to do something with locaiton of legend...
  #  legend("bottomleft", legend = unique.groups[[jj]],
  #        col=cols[[jj]], lty = c(0,0),pch=c(15,19),lwd = c(1.5,1.5),cex = cex_val)
  #}
  
  #par(oldpar)             # return to old par settings
  
  ## Stop saving to current file
  #if (saveplot == TRUE) { 
  #  dev.off()
  #}
  
  ########################################################
  
  ## Not needed as a plot with sample names has the group data as well (by colour)
  ## Plot PCA scores with group names
  ##pic_gen(pca_mat,
  ##   plot_title=if (!is.null(main)) main else "PCA Score Plot\nGroups",
  ##    plot_labels=group_list,
  ##    x_label=paste("PC", x.axis, " (", x_percent, "%)", sep=""),
  ##    y_label=paste("PC", y.axis," (", y_percent, "%)",  sep=""),
  ##    cols_used=cols_used,
  ##    pch_used=pch_used
  ##)
  
  ##check if plot needs to be saved
  #if (saveplot == TRUE) {
  #  if (savetype == "pdf") {
  #    savef(savenames[4],width = 10, height = 9)
  #  }else {
  #    savef(savenames[4],width = 840, height = 720)
  #  }
  #  #wait for key stroke to show next plot
  #} else {
  #  cat("Hit <Return> to see next plot: ")
  #  line <- readline()
  #}
  
  #eigen_mat <- cbind(eigenvecs[, x.axis],eigenvecs[, y.axis])
  ## Loadings plot
  #pic_gen(eigen_mat,
  #        "PCA Loading Plot",
  #        x_label=paste("PC", x.axis," (", x_percent, "%)", sep=""),
  #        y_label=paste("PC", y.axis," (", y_percent, "%)", sep=""),
  #        cols_used="black", cex_val=cex_val
  #        #text_on=FALSE
  #)
  ## Stop saving to current file
  #if (saveplot == TRUE) { 
  #  dev.off()
  #}
}

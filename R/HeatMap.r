#' Heat map
#' 
#' Produces an interactive or a non-interactive heat map of a metabolomics data matrix optionally clustered
#' according to specified methods
#' 
#' 
#' @param featuredata A data frame in the met data format. This should have sample
#' names in the first column to be read as row names and the metabolomics variables 
#' in the remaining columns. 
#' @param groupdata A data frame or a table with optional sample names in the 
#' first column to be read as row names and group names in the following column.
#' @param saveplot A logical indication whether to save the plot produced.
#' @param plotname Name of the output file if the file is to be saved.
#' @param savetype The required format for the plot to be saved in. Threre is a
#' choice of \code{"png","bmp","jpeg","tiff","pdf"} type files.
#' @param interactiveplot A logical indication whether an interactive plot
#' is to be shown.
#' @param saveinteractiveplot A logical indication whether to save the interactive plot as
#' an \code{"html"} file.
#' @param return.interactive A logical indication whether an interactive plot should be returned as 
#' a variable by the function.
#' @param numeric.mets A logical indication whether metabolite names are numeric. If \code{TRUE}, "m"
#' is added before each metabolite name displyed on the graph.
#' @param colramp A vector containing (in order), the desired number of color elements in the panel,
#' color to use for the lowest, color to use for the highest for the non-interactive plot.
#' @param scale A character indicating if the values should be scaled
#' metabolite-wise ("\code{row}") or group-wise ("\code{column}").
#' @param dendrogram A character indicating whether to draw "\code{none}",
#' "\code{row}", "\code{column}" or "\code{both}" dendrograms for non-interactive plots.
#' @param distmethod The distance measure to be used. This must be one of
#' "\code{euclidean}", "\code{maximum}", "\code{manhattan}", "\code{canberra}",
#' "\code{binary}" or "\code{minkowski}".
#' @param aggmethod The agglomeration method to be used. This should be one of
#' "\code{ward}", "\code{single}", "\code{complete}", "\code{average}",
#' "\code{mcquitty}", "\code{median}" or "\code{centroid}".
#' @param margins A numeric vector of length 2 containing the margins for group
#' and metabolite names, respectively.
#' @param key A logical indicating whether a colour key must be drawn.
#' @param keysize A numeric indicating the size of the key.
#' @param cexRow A numeric indicating the size of the metabolite names.
#' @param ColSideColors A character vector indicating the colours different
#' groups.
#' @param ... Arguments to be passed on to other methods.
#' @author Alysha M De Livera, Gavriel Olshansky
#' @seealso \code{\link[graphics]{par}}, \code{heatmap.2}.
#' @examples
#' 
#'     data(mixdata)  #unadjusted data
#'     HeatMap(mixdata$featuredata,mixdata$sampledata[,1], 
#'              saveplot = FALSE, 
#'             interactiveplot = TRUE, scale = "row", 
#'             dendrogram = "none", colramp=c(75,"magenta","green"))
#' 
#' @export HeatMap
HeatMap <- function (featuredata, groupdata, saveplot=FALSE, plotname="heatmap",
                     savetype= c("png","bmp","jpeg","tiff","pdf"),
                     interactiveplot = TRUE, saveinteractiveplot=FALSE,   #GO 
                     return.interactive = FALSE, numeric.mets = FALSE,
                     colramp = c(75, "magenta", "green"), 
                     scale = c("row", "column", "none"),
                     dendrogram = c("column", "row", "both", "none"),
                     distmethod = "euclidean", aggmethod = "complete", 
                     margins = c(5, 5), key = TRUE, keysize = 1.5, cexRow = 0.5, 
                     ColSideColors = NULL, ...) 
{
  
    # Record Par setting and restore them when exiting
    .parold <- par(no.readonly = TRUE)
    on.exit(par(.parold))
    
    #add m to name if needed, converts met name to string
    if(numeric.mets==TRUE){
      colnames(featuredata) <- paste0("m",colnames(featuredata))
    }
    
    # Convert input to inputdata type
    inputdata <- ToInputdata(featuredata,groupdata)
  
    # Create dendogram
    scale <- match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    
    # Get the saving format
    savetype <- match.arg(savetype)
    
    #add m to name if needed, converts met name to string
    if(numeric.mets==TRUE){
      colnames(featuredata) <- paste0("m",colnames(featuredata))
    }
    
    # Create groups information, colour scale
    groups <- factor(inputdata[, 1], levels = unique(inputdata[, 1]))
    unique.groups <- levels(groups)
    if (is.null(ColSideColors)) {
      cols <- ColList(length(unique.groups))
      ColSideColors <- c(rep(NA, length(rownames(inputdata))))
      for (ii in 1:length(inputdata[, 1])) {
        selected <- which(unique.groups == inputdata[, 1][ii])
        ColSideColors[ii] <- cols[selected]
        }
    }
    
    #Added 19092017
    if(dendrogram=="none")
      Rowv<-FALSE
    else
      Rowv<-TRUE
      
    
    # Create heatmap
    inputdata <- t(editcolnames(inputdata)[, -1])
    p <- heatmap_fn(inputdata, saveplot, plotname, savetype, interactiveplot, saveinteractiveplot,
       return.interactive,
       col = colorpanel(n=as.numeric(colramp[1]), low=colramp[2],  high=colramp[3]), 
       scale = scale, dendrogram = dendrogram, 
       margins = margins, key = key, symkey = FALSE, density.info = "none", 
       lvtrace = "none", ColSideColors = ColSideColors, cexRow = cexRow, 
       keysize = keysize, distfun = function(x) dist(x, method = distmethod), 
       hclustfun = function(x) hclust(x, method = aggmethod), Rowv=Rowv,
       ...)
    
    if (return.interactive ==TRUE){
      return(p)
    }
}

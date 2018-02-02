#' Dendrogram
#' 
#' Performs hierarchical cluster analysis given a distance measure 
#' and an agglomeration method, and produces a dendrogram.
#' 
#' 
#' @param featuredata A data frame in the featuredata format. This should have sample
#' names in the first column to be read as row names and the metabolomics variables 
#' in the remaining columns. 
#' @param groupdata A data frame or a table with optional sample names in the 
#' first column to be read as row names and group names in the following column.
#' @param saveplot A logical indication whether to save the dendrogram produced.
#' @param plotname Name of the output file if the file is to be saved.
#' @param savetype The required format for the plot to be saved in. Threre is a
#' choice of \code{"png","bmp","jpeg","tiff","pdf"} type files.
#' @param distmethod The distance measure to be used. This must be one of
#' "\code{euclidean}", "\code{maximum}", "\code{manhattan}", "\code{canberra}",
#' "\code{binary}" or "\code{minkowski}".
#' @param aggmethod The agglomeration method to be used. This should be one of
#' "\code{ward}", "\code{single}", "\code{complete}", "\code{average}",  #ward -> ward.D
#' "\code{mcquitty}", "\code{median}" or "\code{centroid}".
#' @param main Plot title.
#' @param cex A numerical value giving the amount by which plotting text and
#' symbols should be magnified relative to the default.
#' @param clust A logical indicating whether the results from heirarchical clustering should be grouped
#' @param rect A logical indicatng whether rectanges should be drawn highlighting the groups from clust above
#' @param nclust The desired number of clusters for clust or rect. 
#' @param height The desifed height to obtain clusters for clust or rect. Either nclust or height must be supplied for 
#' clust or rect. 
#' @param bordercol If rect=TRUE, a vector with border colors for the rectangles.
#' @param ... Arguments to be passed on to other methods.
#' @return A dendrogram plot and a list containing an object of class `hclust' and a vector with cluster 
#'  membership if clust or rect is set to TRUE.
#' @author Alysha M De Livera, Gavriel Olshansky
#' @seealso \code{\link[stats]{dist}}, \code{\link[stats]{hclust}}.
#' @examples
#' 
#'     data(mixdata) #unadjusted data
#'     Dendrogram(mixdata$featuredata,mixdata$sampledata[,1])
#' 
#' @export Dendrogram
Dendrogram <- function(featuredata, groupdata, saveplot=FALSE, plotname="dendrogram",
                       savetype= c("png","bmp","jpeg","tiff","pdf"),
                       distmethod="manhattan", aggmethod="ward.D", 
                       main="Dendrogram", cex=0.8,
                       clust=FALSE, rect=FALSE,nclust=NULL,height=NULL,bordercol=2, ...)
{
  # Get groups information
  groups <- factor(groupdata, levels=unique(groupdata))
  
  # Rename featuredata
  hcadata <- featuredata
  
  # Prepare distance matrix for samples
  dist_sample <- dist(hcadata, method=distmethod, ...)
  hca_sample <- hclust(dist_sample, method=aggmethod, ...)
  
  smpl_labels<-paste(groups, ':', rownames(hcadata), sep='')
  
  par_defs<-par(mai=par()$mai,
                mar=par()$mar,
                mex=par()$mex,
                mgp=par()$mgp,
                oma=par()$oma,
                omd=par()$omd,
                omi=par()$omi,
                plt=par()$plt,
                usr=par()$usr,
                xaxp=par()$xaxp,
                xpd=par()$xpd,
                yaxp=par()$yaxp
  )
  on.exit(par(par_defs))
  
  lmat<-matrix(c(4,0,0,2,1,3), ncol=2)
  ### These need to be determined properly
  lwid<-c(0.6,4.25)
  lhei<-c(4.25, 0.25, 1)
  col_list<-ColList(length(levels(groups)))
  cols_used<-col_list[groups[hca_sample$order]]
  unique.groups <- levels(groups)
  cols <- ColList(length(unique.groups))
  
  # Show where things are supposed to go
  #dev.new()
  #plotmat<-layout(lmat, lwid, lhei)
  #layout.show(plotmat)
  
  #get format for saving the file if needed
  savetype <- match.arg(savetype)
  
  #save plot if needed
  if (saveplot == TRUE) {
    savef <- match.fun(savetype)
    if (savetype == "pdf") {
      savef(paste(c(plotname,".",savetype),collapse = ""),width = 10, height = 9)
    }else {
      savef(paste(c(plotname,".",savetype),collapse = ""),width = 840, height = 720)
    }
  }
  
  # Do plots
  #dev.new()
  layout(lmat, lwid, lhei)
  # 1) group colours
  # Use rect() for correct placement of group colours (image() width wrong)
  xpos<-seq(0, 1, length.out=length(smpl_labels))
  xwid<-xpos[2]/2
  par(mar=c(0,0,0,1))
  plot(matrix(c(0,0,1,1),nrow=2,byrow=TRUE), type="n", axes=FALSE, 
       xlab="", ylab=""
  )
  for (ii in 1:length(xpos)) {
    rect(xpos-xwid, 0, xpos+xwid, 1, density=NA, col=cols_used)
  }
  
  # 2) Dendro
  par(mar=c(0,0,0,1))
  plot(hca_sample, 
       labels=rep("", length(groups)),        # clear labels here
       hang=(-1),                             # even-length ends
       main=main,                             # no plot title 
       ########### fix (oma? mtext?)
       axes=FALSE,                            # suppress axes (see below)
       xlab="",                               # suppress all axis titles
       ylab="", 
       sub="",
       ...
       
  )
  
  axmap<-par()$usr
  dendyaxp<-par()$yaxp
  # 3) axis labels (x-axis)
  par(cex=cex, mar=c(0,0,0,1), usr=axmap)
  axis(1, 1:length(smpl_labels), labels=reducedlabel(smpl_labels[hca_sample$order],
                                                     50),
       las=2, tick=0, 
       line=2 # if this is not displaced, it is placed _over_ the group colours
  )
  
  # 4) height bar (y-axis)
  par(cex=cex, mar=c(0,1,0,0), usr=axmap, yaxp=dendyaxp)
  ## interestingly, for very large values of height (i.e. 0 to 1.2e+09; from
  ## raw data), this hits a memory limit. 
  ## use function pretty_simp to deal with the problem
  axis(2, (pretty_simp(dendyaxp[1],dendyaxp[2], n=dendyaxp[3]+1)), 
       line=0.5, las=2
  )
  mtext(side = 2, "Height", line = 4)
  
  #Add in the rect
  if(clust | rect){
    if (is.null(nclust) & is.null(height))
      stop("One of nclust or height must be specified")
    clusters<-cutree(tree=hca_sample,k=nclust,h=height)

  } else clusters<-NULL
  if (rect)
    rect.hclust(tree=hca_sample,k=nclust,h=height,border=bordercol,cluster=clusters)
  
  
  #stop saving output
  if (saveplot == TRUE) {
    legend("topright",bg = "white", legend = unique.groups,
           col=cols, lty = c(1,1),lwd = c(2.5,2.5),cex = 2)
    dev.off()
  } else {
    legend("topright", legend = unique.groups,
           col=cols, lty = c(1,1),lwd = c(2.5,2.5),cex = 1.2)
  }
  
  
  par(par_defs)
  
  
  return(list(hca=hca_sample,clusters=clusters))
}


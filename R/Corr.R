#' Computes correlation matrix for a metabolomics dataset or 
#' compares the correlation between two metabolomics datasets
#' 
#' @param featuredata1 A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' @param featuredata2 A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' @param method Must be one of "pearson", "spearman" or "kendall"
#' @param padjmethod p-value adjustment method. Must be one of "holm", 
#' "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" or "none"
#' @param saveoutput A logical indicating whether the results
#' should be saved as a .csv file.
#' @param outputname The name of the output file if the output has to be saved.
#' @return The result is an object of class `results'.
#'  @seealso \code{\link[DiffCorr]{comp.2.cc.fdr}}.
#'  
#' @author Alysha M De Livera, Gavriel Olshansky
#' @references Fukushima, A. Gene (2013) 518, 209-214
#' 
#' @examples
#' data("featuredata_roots")
#' featuredata_roots[featuredata_roots==0]<-NA
#' imp_data<-MissingValues(LogTransform(featuredata_roots)$featuredata)$featuredata
#' Corr( imp_data[c(1:17),], imp_data[c(18:37),])
#' Corr( imp_data[c(1:17),])
#' 
#' @export Corr
Corr <-function(featuredata1=NULL, featuredata2=NULL, 
                method=c("pearson", "kendall", "spearman"),
                padjmethod=c("BH", "holm", "hochberg", "hommel", 
                             "bonferroni", "BY", "fdr", "none"),
                saveoutput=FALSE, 
                outputname=NULL)
{
  # Match the method
  method <- match.arg(method)
  padjmethod <- match.arg(padjmethod)
  
  # If method is not one of the above listed, then stop
  if (!is.element(method, 
                  c("pearson", "kendall", "spearman"))) 
    stop("Invalid method")
  
  
  if (!is.element(padjmethod, 
                  c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                    "fdr", "none")))
    stop("Invalid padjmethod")
  
  if (is.null(featuredata1))
    stop("enter featuredata1")
  alldatacheck(featuredata=featuredata1)
  
  if (is.null(featuredata2)){
    results<-comp.2.cc.fdr_fn2(data=t(featuredata1), 
                               method = method, 
                               p.adjust.methods = padjmethod,
                               threshold = 1) 

  } else{
    alldatacheck(featuredata=featuredata2)
    outdata<-comp.2.cc.fdr_fn( data1=t(featuredata1), 
                               data2=t(featuredata2), 
                               method = method, 
                               p.adjust.methods = padjmethod,
                               threshold = 1) 
    results <- outdata[,c(1,2,8,7,11)]
    
  }
  
  # Generate the output matrix in .csv format
  if (saveoutput) {
    write.csv(results,
              if (!is.null(outputname)) {
                paste(c(outputname, ".csv"), collapse="")
              } else {
                paste(c("CorrResults_", method, ".csv"), collapse="")
              }
    )
  }
  
  output<-list()
  output$results<-results
  return(structure(output, class="results"))
  
}


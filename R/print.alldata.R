print.alldata <- function(x, ...)
{
#    x$featuredata <- editcolnames(x$featuredata)
  if (!is.null(x$featuredata)){
      cat("\nFeature data\n")
    print(dataview(x$featuredata))
  }
  
    if (!is.null(x$sampledata)){
    cat("\nSample data\n")
    print(dataview(x$sampledata))}
    
    if(!is.null(x$metabolitedata)){
    cat("\nMetabolite data\n")
    print(dataview(x$metabolitedata))
    }

    if(!is.null(x$uvdata)){
      cat("\nUV data\n")
      print(dataview(x$uvdata))
    }
    if(!is.null(x$resdata)){
      print(dataview(x$resdata))
    }
    
        }

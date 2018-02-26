#' Edit column names of a metabolomic data matrix
#' 
#' Edits column names of a metabolomic data matrix to remove the letter `X'
#' appearing at the beginning of metabolite names when they begin with a
#' number.
#' 
#' 
#' @param y A data matrix with metabolite names which need to be corrected.
#' @return A data matrix with corrected metabolite names.

editcolnames <- function(y)
{
    colnames(y) <- if (
        length(
            grep("^X[\\d]", colnames(y), perl=TRUE)
        ) != 0
    ) {
        gsub("^X([\\d].*)", "\\1", colnames(y), perl=TRUE)
    } else {
        colnames(y)
    }
    
    return(y)
}

#' Contrast matrix
#' 
#' Generates a contrast matrix with specified contrasts.
#' 
#' 
#' @param contrasts A character vector with specified contrasts.
#' @param levels A character vector or a factor with levels in the design
#' matrix.
#' @return A contrast matrix.
#' @author Alysha M De Livera
#' @seealso \code{\link[limma]{makeContrasts}}.
#' @examples
#' 
#'     ContrastMatrix(contrasts = c("A-B", "B-C"), levels = c("A", "B", "C", "D"))
#' 
#' @export ContrastMatrix
ContrastMatrix<-function(contrasts, levels){

  mat<-makeContrasts(contrasts=contrasts,levels=levels)
  return(mat)
}

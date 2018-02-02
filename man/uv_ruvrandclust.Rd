\docType{data}
\name{uv_ruvrandclust}
\alias{uv_ruvrandclust}
\title{The \code{\link[NormalizeMets:UVdata]{UVdata}} normalized by the RUV random for clustering method}
\description{
  uv_ruvrandclust is a list containing featuredata (normalized data) and uvdata (the removed unwanted variation component). See the NormalizeMets Vignette and De Livera et al 2015.
}
\examples{
  data(uv_ruvrandclust)
  dataview(uv_ruvrandclust$featuredata)
  dataview(uv_ruvrandclust$uvdata)
}
\author{
  Alysha M De Livera, Gavriel Olshansky
}
\keyword{datasets}

\references{
  De Livera, Alysha M De, M. Aho-Sysi, Laurent Jacob, 
  J. Gagnon-Bartch, Sandra Castillo, J.A. Simpson, and Terence P. Speed. 
  2015. Statistical Methods for Handling Unwanted Variation in 
  Metabolomics Data. \emph{Analytical Chemistry} 87 (7). American Chemical Society: 
    3606-3615. 
}

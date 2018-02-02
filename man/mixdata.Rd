\docType{data}
\name{mixdata}
\alias{mixdata}
\title{A designed metabolomics dataset}
\description{
  mixdata consists of featuredata, sampledata and metabolitedata as described by redestig et al 2009.
  A total of 42 samples with dilution mixtures with 44 components, set at three different alternating concentrations, have been run using GC-TOF MS.
}
\examples{
  data(mixdata)
  dataview(mixdata$featuredata)
  dataview(mixdata$metabolitedata)
}
\author{
  Alysha M De Livera, Gavriel Olshansky
}
\keyword{datasets}

\references{
Redestig, H., Fukushima, A., Stenlund, H., Moritz, T., Arita, M., Saito, K.,
Kusano, M. (2009) Compensation for systematic cross-contribution improves
normalization of mass spectrometry based metabolomics data. \emph{Analytical
Chemistry} 81(19): 7974-7980.

}

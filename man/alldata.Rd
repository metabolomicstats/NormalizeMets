\name{alldata-class}
\docType{class}
\alias{alldata-class}
\title{Metabolomic alldata- class}
\description{The container for output from most of the functions in the 
\code{\link{NormalizeMets}} package.}

\details{
  This object class is a list containing the output matrices that are produced by the function. See NormalizeMets Vignette for details.}

\author{Alysha M De Livera, Gavriel Olshansky}

\examples{
     data(alldata_eg)
     dataview(alldata_eg$featuredata)
     dataview(alldata_eg$sampledata)
     dataview(alldata_eg$metabolitedata)
     }
     
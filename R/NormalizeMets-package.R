
#' NormalizeMets package
#' 
#' @description 
#' The NormalizeMets package is a collection of functions for the statistical analysis of
#'  metabolomics data.
#' 
#' Metabolomics data are inevitably subject to a component of unwanted
#' variation, due to factors such as batch effects, matrix effects, and confounding
#' biological variation. This package is a collection of functions designed to implement, 
#' assess, and choose a suitable normalization method for a given metabolomics study.
#' 
#' @details 
#' Run \code{browseVignettes("NormalizeMets")} or 
#' \code{vignette("NormalizeMets_vignette", package = "NormalizeMets")} for more details.
#' 
#' @name NormalizeMets-package
#' @aliases NormalizeMets
#' @docType package
#' @author Alysha M De Livera, Gavriel Olshansky
#' 
#' 
#' Maintainer: Alysha M De Livera <alyshad@unimelb.edu.au>
NULL

if(getRversion() >= "2.15.1")  utils::globalVariables(c("mixdata"))


#' Metabolomics all data - class
#' 
#' Container for outpust produced by the functions in the
#' \code{\link{NormalizeMets}} package.
#' 
#' 
#' @name alldata
#' @docType class
#' @section Slots/Components: This object class is a list structure that has:
#' \describe{ \item{list("featuredata")}{The metabolomics feature matrix.} 
#' \item{list("sampledata")}{The sample specific data.}, 
#' and \item{list("metabolitedata")}{The metabolite specific data.}}
#' @aliases alldata-class
#' @author Alysha M De Livera, Gavriel Olshansky






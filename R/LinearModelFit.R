#' Linear models 
#' 
#' Fit a linear model to each metabolite in a metabolomics data matrix, 
#' and obtain the coefficients, 95% confidence intervals and p-values. 
#' The featuredata must be log transformed, but does not have to be normalised a priori as the
#' LinearModelFit function can be used to fit the ruv2 method to 
#' accommodate the unwanted variation in the model. Either ordinary
#' statistics or empirical Bayes statistics can be obtained.
#' 
#' 
#' @param featuredata featuredata A data frame in the featuredata format. 
#'  This is a dataframe with metabolites in columns and samples in rows.
#' Unique sample names should be provided as row names.
#' @param factormat A design matrix for the linear model, 
#' consisting of biological factors of interest.
#' @param contrastmat An optional contrast matrix indicating which contrasts
#' need to be tested to answer the biological question of interest.
#' @param ruv2 A logical indicating whether to use the \code{ruv2} method for
#' removing unwanted variation.
#' @param k If \code{ruv2=TRUE}, the number of unwanted variation factors to be
#' included in the model.
#' @param qcmets If \code{ruv2=TRUE}, a vector indicating which metabolites should be used as the
#' internal, external standards or other quality control metabolites. 
#' @param moderated A logical indicating whether moderated statistics should be
#' computed.
#' @param padjmethod A character string specifying p value adjustment method
#' for multiple comparisons. Must be one of "\code{bonferroni}", "\code{holm}"
#' (Holm 1979), "\code{hochberg}" (Hochberg 1988), "\code{hommel}" (Hommel
#' 1988), "\code{BH}" (Benjamini and Hochberg 1995), "\code{BY}" (Benjamini and
#' Yekutieli 2001), or "\code{none}". The default method is set to "\code{BH}".
#' @param ci_alpha Significance level for the confidence intervals.
#' @param saveoutput A logical indicating whether the normalised data matrix
#' should be saved as a csv file.
#' @param outputname The name of the output file if the output has to be saved.
#' @param ... further arguments to be passed to or from methods.
#' @return The result is an object of class
#' \code{\link[limma:marraylm]{MArrayLM}}, containing F statistics, t
#' statistics, corresponding confidence intervals, and adjusted and unadjusted p-values 
#' (see De Livera \emph{et al}., 2012a, 2012b). If \code{moderated=TRUE}, moderated statistics 
#' will be computed by empirical Bayes shrinkage of the standard errors towards a common value (Loennstedt
#' \emph{et al} 2002; Smyth 2004).
#' @author Alysha M De Livera, Gavriel Olshansky
#' @seealso \code{\link[limma]{eBayes}}, \code{\link{ContrastMatrix}}
#' @references Benjamini, Y., Hochberg, Y. (1995) Controlling the false
#' discovery rate: a practical and powerful approach to multiple testing.
#' \emph{Journal of the Royal Statistical Society. Series B (Methodological)}
#' 57(1): 289-300.
#' 
#' Benjamini, Y., Yekutieli, D. (2001) The Control of the False Discovery Rate
#' in Multiple Testing under Dependency. \emph{The Annals of Statistics} 29(4):
#' 1165-1188.
#' 
#' De Livera, A. M., Dias, D. A., De Souza, D., Rupasinghe, T., Pyke, J., Tull,
#' D., Roessner, U., McConville, M., Speed, T. P. (2012a) Normalising and
#' integrating metabolomics data. \emph{Analytical Chemistry} 84(24):
#' 10768-10776.
#' 
#' De Livera, Alysha M De, M. Aho-Sysi, Laurent Jacob, 
#' J. Gagnon-Bartch, Sandra Castillo, J.A. Simpson, and Terence P. Speed. 
#' 2015. Statistical Methods for Handling Unwanted Variation in 
#' Metabolomics Data. \emph{Analytical Chemistry} 87 (7). American Chemical Society: 
#' 3606-3615.
#' 
#' De Livera, A.M., Olshansky, M., Speed, T. P. (2013) Statistical analysis of
#' metabolomics data. \emph{Methods in Molecular Biology} (Clifton, N.J.) 1055: 
#' 291-307. 
#' 
#' Gagnon-Bartsch, Johann A., Speed, T. P. (2012) Using control genes to
#' correct for unwanted variation in microarray data. \emph{Biostatistics}
#' 13(3): 539-552.
#' 
#' Hochberg, Y. (1988) A sharper Bonferroni procedure for multiple tests of
#' significance. \emph{Biometrika} 75(4): 800-802.
#' 
#' Holm, S. (1979) A simple sequentially rejective multiple test procedure.
#' \emph{Scandinavian Journal of Statistics} 6(2): 65-70.
#' 
#' Hommel, G. (1988) A stagewise rejective multiple test procedure based on a
#' modified Bonferroni test. \emph{Biometrika} 75(2): 383-386.
#' 
#' Loennstedt, I., Speed, T. P. (2002) Replicated microarray data.
#' \emph{Statistica Sinica} 12: 31-46.
#' 
#' Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing
#' differential expression in microarray experiments. \emph{Statistical
#' Applications in Genetics and Molecular Biology} 3(1): 3.
#' @examples
#' data("alldata_eg")
#' featuredata_eg<-alldata_eg$featuredata
#' dataview(featuredata_eg)
#' sampledata_eg<-alldata_eg$sampledata
#' dataview(sampledata_eg)
#' metabolitedata_eg<-alldata_eg$metabolitedata
#' dataview(metabolitedata_eg)
#'
#' logdata <- LogTransform(featuredata_eg)
#' dataview(logdata$featuredata)
#' imp <-  MissingValues(logdata$featuredata,sampledata_eg,metabolitedata_eg,
#'                      feature.cutof=0.8, sample.cutoff=0.8, method="knn")
#' dataview(imp$featuredata)
#'
#' #Linear model fit using unadjusted data
#' factormat<-model.matrix(~gender +Age +bmi, sampledata_eg)
#' unadjustedFit<-LinearModelFit(featuredata=imp$featuredata,
#'                              factormat=factormat,
#'                              ruv2=FALSE)
#' unadjustedFit
#'
#' #Linear model fit using `is' normalized data 
#' Norm_is <-NormQcmets(imp$featuredata, method = "is", 
#'                     isvec = imp$featuredata[,which(metabolitedata_eg$IS ==1)[1]])
#' isFit<-LinearModelFit(featuredata=Norm_is$featuredata,
#'                      factormat=factormat,
#'                      ruv2=FALSE)
#' isFit
#'
#' #Linear model fit with ruv-2 normalization
#' ruv2Fit<-LinearModelFit(featuredata=imp$featuredata,
#'                        factormat=factormat,
#'                        ruv2=TRUE,k=2,
#'                        qcmets = which(metabolitedata_eg$IS ==1))
#' ruv2Fit
#'
#' #Linear model fit with ruv-2 normalization, obtaining moderated statistics
#' ruv2FitMod<-LinearModelFit(featuredata=imp$featuredata,
#'                           factormat=factormat,
#'                           ruv2=TRUE,k=2,moderated = TRUE,
#'                           qcmets = which(metabolitedata_eg$IS ==1))
#' ruv2FitMod
#' 
#' @export LinearModelFit
LinearModelFit <- function(featuredata, 
                           factormat=NULL,
                           contrastmat=NULL,
                           ruv2=TRUE,
                           k=NULL, qcmets=NULL, 
                           moderated=FALSE,
                           padjmethod="BH",
                           ci_alpha=0.05,
                           saveoutput=FALSE,
                           outputname="results", ...)
{
  featuredata <- as.matrix(featuredata)
  if (mode(featuredata)!="numeric" ) {
    stop("featuredata must be a numerical data matrix")
  }
  Y <- t(featuredata)
  
  if (ruv2){
    # If no k, get them to enter it
    if (is.null(k)) {
      stop("Please enter the number of unwanted variation factors")
    }
    # If there is no qcmets, get them to enter it
    if (is.null(qcmets)) {
      stop(
        paste("Please enter a vector with columns",
              "corresponding to quality control metabolites"
        )
      )
    }
    
    rzY <- t(Y)
    Wruv2 <- svd(rzY[, qcmets] %*% t(rzY[, qcmets]))$u[, 1:k, drop = FALSE]
    colnames(Wruv2) <- paste("k", c(1:k), sep="")
    designmat <- cbind(factormat, Wruv2)
    
    if (!is.null(contrastmat)) {
      contrastmat <- rbind(
        contrastmat, matrix(0, nrow=k, ncol=ncol(contrastmat))
      )
      rownames(contrastmat) <- colnames(designmat)
    }
  } else {
    designmat <- factormat
  }
  
  #Linear model fit
  fit <- lmFit(Y, design=designmat, ...)
  #Contrast fits if required
  if (!is.null(contrastmat)) {
    fit <- contrasts.fit(fit=fit, contrasts=contrastmat, ...)
  }
  
  #eBayes fit
  ebfit <- eBayes(fit, ...)
  ebfit$metabolites <- ebfit$genes
  ebfit$genes <- ebfit$rank <- ebfit$assign <- NULL
  ebfit$qr <- ebfit$qraux <- ebfit$pivot <- ebfit$tol <- NULL
  ebfit$cov.coefficients <- ebfit$pivot <- ebfit$lods <- NULL
  
  #statistics
  beta <- ebfit$coeff
  if (moderated) {
    Fstat <- ebfit$F
    Fpval <- ebfit$F.p.value
    tstat <- ebfit$t
    tpval <- ebfit$p.value
    se <- beta / tstat
    df<- ebfit$df.residual
  } else {
    tstat <- sweep(
      (ebfit$coef / ebfit$stdev.unscaled), 1, ebfit$sigma, "/"
    )
    tpval <- 2 * pt(-abs(tstat), df=ebfit$df.residual)
    ebfit$t <- tstat
    df <- ebfit$df.residual
    fstat <- classifyTestsF(ebfit, df=df, fstat.only=TRUE)
    Fstat <- as.vector(fstat)
    df1 <- attr(fstat, "df1")
    df2 <- attr(fstat, "df2")
    if (df2[1] > 1e+06) {
      Fpval <- pchisq(df1 * Fstat, df1, lower.tail=FALSE)
    } else {
      Fpval<- pf(Fstat, df1, df2, lower.tail=FALSE)
    }
    se <- beta / tstat
    
    ebfit$F <- Fstat
    ebfit$F.p.value <- Fpval
    ebfit$p.value <- tpval
    ebfit$df.total <- ebfit$s2.post <- ebfit$stdev.unscaled <- NULL
    ebfit$var.prior <-ebfit$proportion <- ebfit$s2.prior <- NULL
    ebfit$df.prior <- NULL
    ebfit$std.error <- se
    ebfit$t <- tstat
    ebfit$df <- df
    ebfit$Amean <- NULL
    ebfit$sigma<-NULL
  }
  
  # ruvmat
  if (ruv2) {
    ebfit$uvmat <- Wruv2
  }
  # adjusted p for t p value
  tpadj <- matrix(NA,ncol=ncol(tstat),nrow=nrow(tstat))
  colnames(tpadj) <- colnames(tpval)
  row.names(tpadj) <- row.names(tpval)
  for (j in 1:ncol(tpadj)) {
    tpadj[,j] <- p.adjust(tpval[,j], method=padjmethod, 
                          n=length(na.omit(tpval[,j]))
    )
  }
  ebfit$adj.p.value <- tpadj 
  
  #Compute confidence intervals
  ci1_mets<- beta- qt(1-ci_alpha/2, df)*se
  ci2_mets<- beta+ qt(1-ci_alpha/2, df)*se
  ebfit$lower.ci<-ci1_mets
  ebfit$upper.ci<-ci2_mets
  
  mat <- data.frame(Fstat,
                    Fpval,
                    p.adjust(Fpval, method=padjmethod, n=length(na.omit(Fpval))),
                    beta, tstat, se, tpval, tpadj, ci1_mets, ci2_mets
  )
  
  colnames(mat) <- c("F stat",
                     "F p value",
                     "Adjusted F p value",
                     paste("coeff", colnames(ebfit$coeff)),
                     paste("t stat", colnames(ebfit$t)),
                     paste("std error", colnames(ebfit$coeff)),
                     paste("t p value", colnames(ebfit$t)),
                     paste("Adjusted t p value", colnames(ebfit$t)),
                     paste("Lower CI", colnames(ebfit$t)),
                     paste("Upper CI", colnames(ebfit$t))
  )
  
  if (saveoutput) {
    write.csv(mat,paste(c(outputname,".csv"),collapse=""))
  }
  
  ebfit$residuals<-t(residuals(ebfit,Y))
  return(fit=ebfit)
}

# NormalizeMets

Software for Normalising and analysing metabolomics data.


### About the package:

Metabolomics data are inevitably subject to a component of unwanted variation, due to factors such as batch effects, matrix effects, and confounding biological variation. The R package contains a collection of functions to aid in the statistical analysis of metabolomic data and can be used assess, select and implement statistical methods for normalizing metabolomics data. The interactive excel interface provides an opportunity to use these functions through a user-friendly interface in excel.


This file contains the Microsoft Excel user interface for the publicly available R package NormalizeMets. The installation of the R pacakge is required for the software to be used.

#### 1. Steps for installing the R package

  #(1) (i)  Install R from the CRAN website here https://cran.r-project.org/bin/windows/base/   
  #(ii) It is also recommended to install RStudio from here https://www.rstudio.com/products/RStudio/#Desktop
	  # In R, '#' refers to a comment. Text on the right of # is ignored by R (the compiler).

  #(2) In Rstudio or Rgui, Install the packages using the code below:
  
  ```
  install.packages( c("GGally", "DiffCorr","plotly", "ggplot2", "htmlwidgets", "MetNorm","AUC", "metabolomics","crmn", "gplots")
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("limma", "Biobase", "pcaMethods", "multtest",  "statTarget", "e1071", "impute"))
  ```


  #installing/loading the package installr:
  ```
  if(!require(installr)) { install.packages("installr"); require(installr)}        
  ```
  
  #Installing pandoc and RGtk2:
  ```
  install.pandoc()
  ```
  ```
  install.packages("RGtk2", depen=T)
  library(RGtk2) 
  ```
  #Click "Install GTK+" when prompted


  #Restart R/Rstudio and confirm GTK is now up and running by reloading the package: 
  ```
  library(RGtk2)
  ```
  
  #(3) Install NormalizeMets package using:
  ```
  install.packages("NormalizeMets", type = "source")
  ```
  #(4)  Load the package using
  ```
  library(NormalizeMets)        # tests if package loads
  ```

#### 2. To use the interactive excel interface, download and open the file ExNormalizeMets.xlsm


#### 3. The vignette and workflow are available on this website as a zip file above.

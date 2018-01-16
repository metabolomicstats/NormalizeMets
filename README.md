# NormalizeMets

Authors: Alysha M De Livera, Gavriel Olshansky <br />
License: GPL (>=2) <br />
Contact: <br />
* Excel: g.olshansky@student.unimelb.edu.au
* R    : alyshad@unimelb.edu.au
<br />
## Software for Normalising and Analysing metabolomics data.


### About the package:

Metabolomics data are inevitably subject to a component of unwanted variation, due to factors such as batch effects, matrix effects, and confounding biological variation. The NormalizeMEts R package contains a collection of functions to aid in the statistical analysis of metabolomic data and can be used assess, select and implement statistical methods for normalizing metabolomics data. The interactive excel interface provides an opportunity to use these functions through a user-friendly interface in excel.

#### Requirements
The software requires the installation of R software (version 3.4.3 or higher). Users who wish to use the Excel interface are required to additionally install Microsoft Excel (2016 or later).

#### Input data format

The input data format consists of three parts: (i) "featuredata" which is the metabolomics data matrix containing all metabolite peak intensities (or concentrations). Unique sample names must be provided as row names and unique metabolite names as column names, (ii) "metabolitedata" contains metabolite-specific information in a separate dataframe. These information can include, but is not limited to, designation of metabolites as internal/external standards, or positive/negative controls. Metabolite names need to be provided as row names, and (iii) "sampledata" is a dataframe that contains sample-specific information. These information can include sample type, order of analysis, factors of interest and other sample-specific data relevant to the analysis. Unique sample names need to be provided as row names.


#### The Excel interface

This file ExNormalizeMets.xlsm contains the Microsoft Excel user interface for the publicly available R package NormalizeMets. The installation of the R pacakge is required for the software to be used.

#### Steps for installing the R package (for use in excel):

#### To view a detailed installation and usage guide, open ExNormalizeMets_manual.pdf

#### To use the interactive excel interface, download and open the file ExNormalizeMets.xlsm

#### The vignette, examples and the workflow are available on this website as a html file in the folder NormalizeMets_vignette.zip.

# Trophoblast

Code for analysis of data for the paper "Human and marmoset trophoblast stem cells differ in signalling requirements and recapitulate divergent modes of trophoblast invasion" by Siriwardena et al. 

## Processing scripts

Processing scripts for datasets can be found in the Scripts folder

## Stand alone scripts

Stand alone scripts for individual analyses can be found in the Rscripts folder

## Spatial modelling

Code for spatial idetity mapping can be found in the SpatialMapping folder.

Usage:

Code for performing spatial projection by correlation can be found in the SpatialMapping folder and requires 3 steps: 
1) After downloading the package additional data must be downloaded from here: 
2) The R script genFinalMappingCorrelations.R will load the count data of the in vitro samples and in vivo referece and generate interemdiate files. Note for ease of use these files have bee uploaded into the Output folder in this repo.
3) The ProjectMarmosetDataCorrelation.m code uses GP models to infer spatial identity mapping of individual cells. Note: this code makes use of the Gaussian process (gpml package) of Rasmussenn and Nickisch (http://gaussianprocess.org/gpml/code/matlab/doc). 


# Trophoblast

Code for analysis of data for the paper "Human and marmoset trophoblast stem cells differ in signalling requirements and recapitulate divergent modes of trophoblast invasion ![image](https://user-images.githubusercontent.com/59876617/224363561-1dad3f19-3f75-4a4e-ac62-78c0a16dae7a.png)
" by Siriwardena et al. 

Usage:

Code for performing spatial projection by correlation can be found in the SpatialMapping folder and requires 3 steps: 
1) After downloading the package additional data must be downloaded from here: 
2) The R script genFinalMappingCorrelations.R will load the count data of the in vitro samples and in vivo referece and generate interemdiate files. Note for ease of use these files have bee uploaded into the Output folder in this repo.
3) The 
This code makes use of the Gaussian process (gpml package) of Rasmussenn and Nickisch (http://gaussianprocess.org/gpml/code/matlab/doc). 


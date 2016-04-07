###########################################################################
############# Script for Bayes Factor test in R #################
###########################################################################

### Clear the workspace
rm(list=ls())

## Read in the lstat files
unconstrained.means<-read.table("Data/MrBayes_Output/parameter_files/unconstrained_parameter_files/MrBayes_uncalibrated_Thalattosuchia_March2016.nex.lstat", header=TRUE)
constrained201.means<-read.table("Data/MrBayes_Output/parameter_files/constrained_201Ma_parameter_files/MrBayes_calibrated_Thalattosuchia_201Ma_March2016.nex.lstat", header=TRUE)
constrained215.means<-read.table("Data/MrBayes_Output/parameter_files/constrained_215Ma_parameter_files/MrBayes_calibrated_Thalattosuchia_March2016.nex.lstat", header=TRUE)

uncon.vs.201<-(-constrained201.means$harmonic_mean[11])/(-unconstrained.means$harmonic_mean[11])
uncon.vs.215<-(-constrained215.means$harmonic_mean[11])/(-unconstrained.means$harmonic_mean[11])
con201.vs.215<-(-constrained201.means$harmonic_mean[11])/(-constrained215.means$harmonic_mean[11])

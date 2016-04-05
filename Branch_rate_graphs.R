###########################################################################
############# Script for plotting rates of morphological change over time #################
###########################################################################

## @knitr constrained_201_treeprocess

### Clear the workspace
rm(list=ls())

### Read in the relative branch rates
constrained201.rates<-read.table("Data/MrBayes_Output/parameter_files/constrained_201Ma_parameter_files/constrained_201Ma_IGRbranchrates.txt", sep="\t", header = TRUE)

### calculate absolute rates based on the median relative rate
constrained201.rates$absolute<-constrained201.rates$Median*1.981174473791735e-03

### calculate the median % change per Ma for each branch
constrained201.rates$percentage<-constrained201.rates$absolute*100

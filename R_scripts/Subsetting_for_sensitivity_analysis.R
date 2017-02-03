###########################################################################
############# Script for Bayes Factor test in R #################
###########################################################################

### Clear the workspace
rm(list=ls())

# Read in the data matrix for subsetting
characters<-read.csv("data/sensitivity_testing/Data_completeness.csv", stringsAsFactors = FALSE)
ages<-read.csv("data/sensitivity_testing/Taxon_ages.csv", stringsAsFactors = FALSE)

for(i in c(10,20,30,40,50)){
  trimmed<-subset(characters,Percentage_missing>i, select=c("Species","Characters"))
  write.table(trimmed, file=paste("data/sensitivity_testing/Subset_characters_",i,".txt",sep=""), sep=" ", col.names = FALSE, row.names = FALSE)
}

for(i in c(10,20,30,40,50)){
  trimmed_characters<-subset(characters,Percentage_missing>i, select=c("Species","Characters"))
  trimmed_ages<-ages[match(trimmed_characters$Species,ages$Taxa),]
  write.csv(trimmed_ages, file=paste("data/sensitivity_testing/Subset_ages_",i,".csv",sep=""), row.names = FALSE)
}
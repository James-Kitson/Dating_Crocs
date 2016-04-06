########################################################################################################
######################## Script for plotting Cratopus trees ############################
########################################################################################################

### Clear the workspace
rm(list=ls())

### open APE
library(ape)
library(phytools)
library(plyr)

### read in the trees
my.trees201<-read.nexus("Data/Consensus_trees/constrained_201Ma_Thallatosuchia_MCC_tree")
my.trees215<-read.nexus("Data/Consensus_trees/constrained_215Ma_Thallatosuchia_MCC_tree")
my.trees.unconstrained<-read.nexus("Data/Consensus_trees/unconstrained_Thallatosuchia_MCC_tree")

my.tree201<-my.trees201
my.tree215<-my.trees215
my.tree.unconstrained<-my.trees.unconstrained


### read in the list of names
name<-read.csv("Data/Metadata/March_2016_Taxon_names.csv")
### read.csv turns text into factors, this gets messy later when plotting
### so make it character data
name<-data.frame(lapply(name, as.character), stringsAsFactors=FALSE)

### the next line uses match to perform the same function as vlookup in excel
my.tree201$tip.label <- (name$Real.Names[match(my.tree201$tip.label,name$Taxa)])
my.tree215$tip.label <- (name$Real.Names[match(my.tree215$tip.label,name$Taxa)])
my.tree.unconstrained$tip.label <- (name$Real.Names[match(my.tree.unconstrained$tip.label,name$Taxa)])

### use cophylo to rotate the MCC tree relative to the consensus tree
obj<-cophylo(my.tree201,my.tree215)
obj2<-cophylo(my.tree.unconstrained,my.tree201)

### plot and look at the tip associations
pdf("Diagrams/branching_comparison_201Mavs215Ma_MCC.pdf",30,20)
plot(obj, cex=0.5)
dev.off()

### plot and look at the tip associations
pdf("Diagrams/branching_comparison_unconstrainedvs201Ma_MCC.pdf",30,20)
plot(obj2, cex=0.5)
dev.off()


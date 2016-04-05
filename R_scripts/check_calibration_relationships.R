########################################################################################################
######################## Script for plotting Cratopus trees ############################
########################################################################################################

## @knitr MLtreeprocess

### Clear the workspace
rm(list=ls())

### open APE
library(ape)
library(phytools)
library(plyr)

### read in the trees
my.trees201<-read.nexus("Data/Consensus_trees/MrBayes_constrained_Thalattosuchia_201Ma_March2016.nex.con.tre")
my.trees215<-read.nexus("Data/Consensus_trees/MrBayes_constrained_Thalattosuchia_215Ma_March2016.nex.con.tre")
my.trees.unconstrained<-read.nexus("Data/Consensus_trees/MrBayes_unconstrained_Thalattosuchia_March2016.nex.con.tre")

my.tree201<-my.trees201[[1]]
my.tree215<-my.trees215[[1]]
my.tree.unconstrained<-my.trees.unconstrained[[1]]


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
pdf("Diagrams/branching_comparison_201Mavs215Ma.pdf",30,20)
plot(obj, cex=0.5)
dev.off()

### plot and look at the tip associations
pdf("Diagrams/branching_comparison_unconstrainedvs201Ma.pdf",30,20)
plot(obj2, cex=0.5)
dev.off()


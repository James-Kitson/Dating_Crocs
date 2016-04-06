###########################################################################
############# Script for plotting rates of morphological change over time #################
###########################################################################

## @knitr constrained_201_treeprocess

### Clear the workspace
rm(list=ls())

### load the ape and phytools libraries
library(ape)
library(phytools)
library(Claddis)

### Read in the relative branch rates
constrained201.rates<-read.table("Data/MrBayes_Output/parameter_files/constrained_201Ma_parameter_files/constrained_201Ma_IGRbranchrates.txt", sep="\t", header = TRUE)

### calculate absolute rates based on the median relative rate
constrained201.rates$absolute<-constrained201.rates$Median*1.981174473791735e-03

### calculate the median % change per Ma for each branch
constrained201.rates$percentage<-constrained201.rates$absolute*100

### Load the tree from the working directory
Bayes.tree<-read.nexus("Data/Consensus_trees/MrBayes_constrained_Thalattosuchia_201Ma_March2016.nex.con.tre")

### conformat=simple in MrBayes produces two trees, one with branch lengths
### and one with node probabilities, all the modification needs to happen
### on the first one so make than an object in its own right
My.tree<-Bayes.tree[[1]]
My.tree<-ladderize(My.tree)
My.tree<-root(My.tree,1)

### scale the edge lengths to the root age of the tree
My.tree$edge.length<-My.tree$edge.length/(max(nodeHeights(My.tree)/227))

### set a root age in the tree object for axisPhylo.
My.tree$root.time<-227

GetNodeAges(My.tree)

plot(My.tree)
axisPhylo()


geoscalePhylo(ladderize(My.tree), cex.age=0.6, cex.ts=0.8, cex.tip=1)

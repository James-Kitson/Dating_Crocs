###########################################################################
############# Script for renaming and plotting trees in R #################
########## Requires MrBayes to use conformat=simple and a csv of ##########
############## the OTU codes and corresponding real names #################
###########################################################################

## @knitr lttcomparisonprocess

### Clear the workspace
rm(list=ls())

### load the ape and phytools libraries
library(ape)
library(phytools)
library(strap)

###########################################################################
######################## Get the data #####################################
###########################################################################

### Load the trees from the working directory
Bayes.tree.constrained201<-read.nexus("Data/Consensus_trees/MrBayes_constrained_Thalattosuchia_201Ma_May2016.nex.con.tre")
Bayes.tree.constrained215<-read.nexus("Data/Consensus_trees/MrBayes_constrained_Thalattosuchia_215Ma_May2016.nex.con.tre")
Bayes.tree.unconstrained<-read.nexus("Data/Consensus_trees/MrBayes_unconstrained_Thalattosuchia_May2016.nex.con.tre")

### Process the 201Ma constrained trees
tree.constrained201<-Bayes.tree.constrained201[[1]]
tree.constrained201<-ladderize(tree.constrained201)
tree.constrained201<-root(tree.constrained201,1)
### scale the edge lengths to the root age of the tree
tree.constrained201$edge.length<-tree.constrained201$edge.length/(max(nodeHeights(tree.constrained201)/227))
### set a root age in the tree object for axisPhylo.
tree.constrained201$root.time<-227
### calculate the values for the lineage through time plot
tree.constrained201_coords<-as.data.frame(ltt.plot.coords(tree.constrained201, backward=TRUE))
### invert the ages as the geoscalePlot plots the stratigraphic axis in the opposite direction to geoscalePhylo
tree.constrained201_coords$time<-tree.constrained201_coords$time*(-1)

### Process the 215Ma constrained trees
tree.constrained215<-Bayes.tree.constrained215[[1]]
tree.constrained215<-ladderize(tree.constrained215)
tree.constrained215<-root(tree.constrained215,1)
### scale the edge lengths to the root age of the tree
tree.constrained215$edge.length<-tree.constrained215$edge.length/(max(nodeHeights(tree.constrained215)/227))
### set a root age in the tree object for axisPhylo.
tree.constrained215$root.time<-227
### calculate the values for the lineage through time plot
tree.constrained215_coords<-as.data.frame(ltt.plot.coords(tree.constrained215, backward=TRUE))
### invert the ages as the geoscalePlot plots the stratigraphic axis in the opposite direction to geoscalePhylo
tree.constrained215_coords$time<-tree.constrained215_coords$time*(-1)

### Process the unconstrained trees
tree.unconstrained<-Bayes.tree.unconstrained[[1]]
tree.unconstrained<-ladderize(tree.unconstrained)
tree.unconstrained<-root(tree.unconstrained,1)
### scale the edge lengths to the root age of the tree
tree.unconstrained$edge.length<-tree.unconstrained$edge.length/(max(nodeHeights(tree.unconstrained)/227))
### set a root age in the tree object for axisPhylo.
tree.unconstrained$root.time<-227
### calculate the values for the lineage through time plot
tree.unconstrained_coords<-as.data.frame(ltt.plot.coords(tree.unconstrained, backward=TRUE))
### invert the ages as the geoscalePlot plots the stratigraphic axis in the opposite direction to geoscalePhylo
tree.unconstrained_coords$time<-tree.unconstrained_coords$time*(-1)

## @knitr lttcomparisonplot

###################################
### THIS CAN TOTALLY BE DONE BETTER USING geoscalePlot AND ltt.plot.coords
##################################

geoscalePlot(tree.unconstrained_coords$time,tree.unconstrained_coords$N,
             units=c("Period", "Epoch"),
             boxes= "Epoch", type="l",
             data.lim = c(0,70),
             label= "Number of extant lineages")
lines(tree.constrained215_coords$time,tree.constrained215_coords$N, lty=2)
lines(tree.constrained201_coords$time,tree.constrained201_coords$N, lty=3)


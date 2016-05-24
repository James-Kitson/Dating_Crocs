###########################################################################
############# Script for renaming and plotting trees in R #################
########## Requires MrBayes to use conformat=simple and a csv of ##########
############## the OTU codes and corresponding real names #################
###########################################################################

## @knitr constrained_201_treeprocess

### Clear the workspace
rm(list=ls())

### load the ape and phytools libraries
library(ape)
library(phytools)
library(strap)

###########################################################################
######################## Get the data #####################################
###########################################################################


### Load the tree from the working directory
Bayes.tree<-read.nexus("Data/Consensus_trees/MrBayes_constrained_Thalattosuchia_201Ma_May2016.nex.con.tre")
### conformat=simple in MrBayes produces two trees, one with branch lengths
### and one with node probabilities, all the modification needs to happen
### on the first one so make than an object in its own right
tree<-Bayes.tree[[1]]
tree<-ladderize(tree)
tree<-root(tree,1)
### read in the list of names
name<-read.csv("Data/Metadata/March_2016_Taxon_names.csv")
### read.csv turns text into factors, this gets messy later when plotting
### so make it character data
name<-data.frame(lapply(name, as.character), stringsAsFactors=FALSE)
########################################################

### replace the codes with informative names
tree.rename<-tree
### the next line uses match to perform the same function as vlookup in excel
tree.rename$tip.label <- (name$Real.Names[match(tree.rename$tip.label,name$Taxa)])

### Convert node lables to numeric and round to two dp
tree.rename$node.label<-as.numeric(tree.rename$node.label)
tree.rename$node.label<-round(tree.rename$node.label,digits=2)

### scale the edge lengths to the root age of the tree
tree.rename$edge.length<-tree.rename$edge.length/(max(nodeHeights(tree.rename)/227))
### set a root age in the tree object for axisPhylo.
tree.rename$root.time<-227

## @knitr constrained_201_treeplot

#########################################################
### Plot the main combined tree with support values
#pdf(file=paste(out,file.names[[x]],"_",title,".pdf",sep=""), 30, 30)

geoscalePhylo(ladderize(tree.rename), cex.age=0.6, cex.ts=0.8, cex.tip=0.6, quat.rm=TRUE, units=c("Period", "Epoch"), boxes= "Epoch")

#dev.off()
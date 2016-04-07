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

###########################################################################
######################## Get the data #####################################
###########################################################################

### Load the trees from the working directory
Bayes.tree.constrained201<-read.nexus("Data/Consensus_trees/MrBayes_constrained_Thalattosuchia_201Ma_March2016.nex.con.tre")
Bayes.tree.constrained215<-read.nexus("Data/Consensus_trees/MrBayes_constrained_Thalattosuchia_215Ma_March2016.nex.con.tre")
Bayes.tree.unconstrained<-read.nexus("Data/Consensus_trees/MrBayes_unconstrained_Thalattosuchia_March2016.nex.con.tre")

### Process the 201Ma constrained trees
tree.constrained201<-Bayes.tree.constrained201[[1]]
tree.constrained201<-ladderize(tree.constrained201)
tree.constrained201<-root(tree.constrained201,1)
### scale the edge lengths to the root age of the tree
tree.constrained201$edge.length<-tree.constrained201$edge.length/(max(nodeHeights(tree.constrained201)/227))
### set a root age in the tree object for axisPhylo.
tree.constrained201$root.time<-227

### Process the 215Ma constrained trees
tree.constrained215<-Bayes.tree.constrained215[[1]]
tree.constrained215<-ladderize(tree.constrained215)
tree.constrained215<-root(tree.constrained215,1)
### scale the edge lengths to the root age of the tree
tree.constrained215$edge.length<-tree.constrained215$edge.length/(max(nodeHeights(tree.constrained215)/227))
### set a root age in the tree object for axisPhylo.
tree.constrained215$root.time<-227

### Process the unconstrained trees
tree.unconstrained<-Bayes.tree.unconstrained[[1]]
tree.unconstrained<-ladderize(tree.unconstrained)
tree.unconstrained<-root(tree.unconstrained,1)
### scale the edge lengths to the root age of the tree
tree.unconstrained$edge.length<-tree.unconstrained$edge.length/(max(nodeHeights(tree.unconstrained)/227))
### set a root age in the tree object for axisPhylo.
tree.unconstrained$root.time<-227

### Read in the geological epoch data
epochs<-read.csv("Data/Metadata/epochs_and_colours.csv", stringsAsFactors = F)

### set up a vector of the epoch colours based on Commission for the Geological Map of the World guidelines
legend.cols<-rgb(red=epochs$Red, green=epochs$Green, blue=epochs$Blue, alpha=125, maxColorValue = 255)

## @knitr lttcomparisonplot

###################################
### THIS CAN TOTALLY BE DONE BETTER USING geoscalePlot AND ltt.plot.coords
##################################


### plot the unconstrained ltt
ltt.plot(tree.unconstrained, xlab="Time (Ma)", ylab="Extant lineages", ylim=c(0,70), lty=1)
ltt.lines(tree.constrained215, lty=2)
ltt.lines(tree.constrained201, lty=3)
title(main="Thallatosuchian root unconstrained")
### plot the epochs
for(i in 1:length(epochs$Starting)){
  polygon(x=c(-epochs$Starting[i],-epochs$Starting[i],-epochs$Ending[i],-epochs$Ending[i]),
          y = c(0,70,70,0), border=NA, col =legend.cols[i])
}
unconsltt<-as.data.frame(ltt.plot.coords(tree.unconstrained,backward = FALSE))
geoscalePlot(ages=unconsltt$time, data=unconsltt$N)

legend("topright", title="Epoch", inset=0.005, legend = epochs$Stage,
       fill =legend.cols,
       cex=0.5,
       bg = "white")
legend("right",legend=c("Unconstrained Thallatosuchian\nroot","Thallatosuchian root\n205Ma-215Ma","Thallatosuchian root\n19Ma-201Ma"),
       lty=c(1,2,3),
       inset=0.005,
       cex=0.5,
       bg = "white",
       y.intersp = 2)
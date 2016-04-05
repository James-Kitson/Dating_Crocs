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

### Process the 215Ma constrained trees
tree.constrained215<-Bayes.tree.constrained215[[1]]
tree.constrained215<-ladderize(tree.constrained215)
tree.constrained215<-root(tree.constrained215,1)

### Process the unconstrained trees
tree.unconstrained<-Bayes.tree.unconstrained[[1]]
tree.unconstrained<-ladderize(tree.unconstrained)
tree.unconstrained<-root(tree.unconstrained,1)

### make an offset for the axis as R won't draw it from the tip to the root. The offset is a negative starting point for the axis equivalent to the
### rounding up we do at the root end of the axis i.e. if we round 227 Mya to 230 Mya then we need to offset by minus 3Ma of distance measured in
### branch lengths. To do this we divide the root height by the root age and multiply by -3.
offset.unconstrained<-3*(max(nodeHeights(tree.unconstrained)/227))
offset.constrained201<-3*(max(nodeHeights(tree.constrained201)/227))
offset.constrained215<-3*(max(nodeHeights(tree.constrained215)/227))

### Read in the geological epoch data
epochs<-read.csv("Data/Metadata/epochs_and_colours.csv", stringsAsFactors = F)

### calculate the x co-ordinates for the constrained tree polygons on the 201Ma constrained tree
epochs$constrained201.start<-(-max(nodeHeights(tree.constrained201)+offset.constrained201)/23)*(epochs$Starting/10)
epochs$constrained201.end<-(-max(nodeHeights(tree.constrained201)+offset.constrained201)/23)*(epochs$Ending/10)

### calculate the x co-ordinates for the constrained tree polygons on the 201Ma constrained tree
epochs$constrained215.start<-(-max(nodeHeights(tree.constrained215)+offset.constrained215)/23)*(epochs$Starting/10)
epochs$constrained215.end<-(-max(nodeHeights(tree.constrained215)+offset.constrained215)/23)*(epochs$Ending/10)

### calculate the x co-ordinates for the unconstrained tree polygons
epochs$unconstrained.start<-(-max(nodeHeights(tree.unconstrained)+offset.unconstrained)/23)*(epochs$Starting/10)
epochs$unconstrained.end<-(-max(nodeHeights(tree.unconstrained)+offset.unconstrained)/23)*(epochs$Ending/10)

### set up a vector of the epoch colours based on Commission for the Geological Map of the World guidelines
legend.cols<-rgb(red=epochs$Red, green=epochs$Green, blue=epochs$Blue, alpha=125, maxColorValue = 255)

## @knitr lttcomparisonplot

par(mfrow=c(3,1),
    mai=c(1, 1, 0.2, 1))

### plot the unconstrained ltt
ltt.plot(tree.unconstrained, xaxt="n", xlab="Time (Ma)", ylab="Extant lineages", ylim=c(0,70))
title(main="Thallatosuchian root unconstrained")
### plot the epochs
for(i in 1:length(epochs$Starting)){
  polygon(x=c(epochs$unconstrained.start[i],epochs$unconstrained.start[i],epochs$unconstrained.end[i],epochs$unconstrained.end[i]),
          y = c(0,70,70,0), border=NA, col =legend.cols[i])
}
axis(side=1,cex.axis=0.8, padj=1,at=seq(from=-(max(nodeHeights(tree.unconstrained)+offset.unconstrained)), to=0, by=(max(nodeHeights(tree.unconstrained))+offset.unconstrained)/23), labels=seq(230,0,by=-10))
legend("topright", title="Epoch", inset=0.005, legend = epochs$Stage,
       fill =legend.cols,
       cex=0.5,
       bg = "white")

### plot the ltt for the analysis constrained at 201Ma
ltt.plot(tree.constrained201, xaxt="n", xlab="Time (Ma)", ylab="Extant lineages", ylim=c(0,70))
title(main="Thallatosuchian root constrained between 191Ma and 201Ma")
for(i in 1:length(epochs$Starting)){
  polygon(x=c(epochs$constrained201.start[i],epochs$constrained201.start[i],epochs$constrained201.end[i],epochs$constrained201.end[i]),
          y = c(0,70,70,0), border=NA, col =legend.cols[i])
}
axis(side=1,cex.axis=0.8, padj=1,at=seq(from=-(max(nodeHeights(tree.constrained201)+offset.constrained201)), to=0, by=(max(nodeHeights(tree.constrained201))+offset.constrained201)/23), labels=seq(230,0,by=-10))

### plot the ltt for the analysis constrained at 215Ma
ltt.plot(tree.constrained215, xaxt="n", xlab="Time (Ma)", ylab="Extant lineages", ylim=c(0,70))
title(main="Thallatosuchian root constrained between 205Ma and 215Ma")
for(i in 1:length(epochs$Starting)){
  polygon(x=c(epochs$constrained215.start[i],epochs$constrained215.start[i],epochs$constrained215.end[i],epochs$constrained215.end[i]),
          y = c(0,70,70,0), border=NA, col =legend.cols[i])
}
axis(side=1,cex.axis=0.8, padj=1,at=seq(from=-(max(nodeHeights(tree.constrained215)+offset.constrained215)), to=0, by=(max(nodeHeights(tree.constrained215))+offset.constrained215)/23), labels=seq(230,0,by=-10))

#dev.off()

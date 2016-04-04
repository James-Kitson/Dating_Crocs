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
Bayes.tree.calibrated<-read.nexus("Data/Consensus_trees/MrBayes_calibrated_Thalattosuchia_March2016.nex.con.tre")
Bayes.tree.uncalibrated<-read.nexus("Data/Consensus_trees/MrBayes_uncalibrated_Thalattosuchia_March2016.nex.con.tre")

### Process the calibrated trees
tree.calibrated<-Bayes.tree.calibrated[[1]]
tree.calibrated<-ladderize(tree.calibrated)
tree.calibrated<-root(tree.calibrated,1)

### Process the uncalibrated trees
tree.uncalibrated<-Bayes.tree.uncalibrated[[1]]
tree.uncalibrated<-ladderize(tree.uncalibrated)
tree.uncalibrated<-root(tree.uncalibrated,1)

### make an offset for the axis as R won't draw it from the tip to the root. The offset is a negative starting point for the axis equivalent to the
### rounding up we do at the root end of the axis i.e. if we round 227 Mya to 230 Mya then we need to offset by minus 3Ma of distance measured in
### branch lengths. To do this we divide the root height by the root age and multiply by -3.
offset.uncalibrated<-3*(max(nodeHeights(tree.uncalibrated)/227))
offset.calibrated<-3*(max(nodeHeights(tree.calibrated)/227))

### Read in the geological epoch data
epochs<-read.csv("Data/Metadata/epochs_and_colours.csv", stringsAsFactors = F)

### calculate the x co-ordinates for the calibrated tree polygons
epochs$calibrated.start<-(-max(nodeHeights(tree.calibrated)+offset.calibrated)/23)*(epochs$Starting/10)
epochs$calibrated.end<-(-max(nodeHeights(tree.calibrated)+offset.calibrated)/23)*(epochs$Ending/10)

### calculate the x co-ordinates for the uncalibrated tree polygons
epochs$uncalibrated.start<-(-max(nodeHeights(tree.uncalibrated)+offset.uncalibrated)/23)*(epochs$Starting/10)
epochs$uncalibrated.end<-(-max(nodeHeights(tree.uncalibrated)+offset.uncalibrated)/23)*(epochs$Ending/10)

### set up a vector of the epoch colours based on Commission for the Geological Map of the World guidelines
legend.cols<-rgb(red=epochs$Red, green=epochs$Green, blue=epochs$Blue, alpha=125, maxColorValue = 255)

## @knitr lttcomparisonplot

par(mfrow=c(2,1),
    mai=c(1, 1, 0.2, 1))
ltt.plot(tree.uncalibrated, xaxt="n", xlab="Time (Ma)", ylab="Extant lineages", ylim=c(0,70))
title(main="Thallatosuchian root unconstrained")
### plot the epochs
for(i in 1:length(epochs$Starting)){
  polygon(x=c(epochs$uncalibrated.start[i],epochs$uncalibrated.start[i],epochs$uncalibrated.end[i],epochs$uncalibrated.end[i]),
          y = c(0,65,65,0), border=NA, col =legend.cols[i])
}
axis(side=1,cex.axis=0.8, padj=1,at=seq(from=-(max(nodeHeights(tree.uncalibrated)+offset.uncalibrated)), to=0, by=(max(nodeHeights(tree.uncalibrated))+offset.uncalibrated)/23), labels=seq(230,0,by=-10))
legend("topright", title="Epoch", inset=0.005, legend = epochs$Stage,
       fill =legend.cols,
       cex=0.5,
       bg = "white")

ltt.plot(tree.calibrated, xaxt="n", xlab="Time (Ma)", ylab="Extant lineages", ylim=c(0,70))
title(main="Thallatosuchian root constrained between 205Ma and 215Ma")
for(i in 1:length(epochs$Starting)){
  polygon(x=c(epochs$calibrated.start[i],epochs$calibrated.start[i],epochs$calibrated.end[i],epochs$calibrated.end[i]),
          y = c(0,65,65,0), border=NA, col =legend.cols[i])
}
axis(side=1,cex.axis=0.8, padj=1,at=seq(from=-(max(nodeHeights(tree.calibrated)+offset.calibrated)), to=0, by=(max(nodeHeights(tree.calibrated))+offset.calibrated)/23), labels=seq(230,0,by=-10))
#dev.off()

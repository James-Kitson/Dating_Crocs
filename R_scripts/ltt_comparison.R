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


### Load the tree from the working directory
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

epochs<-read.csv("Data/Metadata/epochs_and_colours.csv")
epochs$calibrated.start<-(-max(nodeHeights(tree.calibrated)+offset.calibrated)/23)*(epochs$Starting/10)
epochs$calibrated.end<-(-max(nodeHeights(tree.calibrated)+offset.calibrated)/23)*(epochs$Ending/10)

## @knitr lttcomparisonplot

par(mfrow=c(2,1))
ltt.plot(tree.uncalibrated, xaxt="n", xlab="Time (Ma)", ylab="Extant lineages")
legend(x=-0.08, y=60, legend ="uncalibrated Thallatosuchia",col ="black", lty=1, bty="n", y.intersp = 1.1)
abline(v=c(-0.3,-0.2))
axis(side=1,cex.axis=1.0, padj=1,at=seq(from=-(max(nodeHeights(tree.uncalibrated)+offset.uncalibrated)), to=0, by=(max(nodeHeights(tree.uncalibrated))+offset.uncalibrated)/23), labels=seq(230,0,by=-10))
ltt.plot(tree.calibrated,xaxt="n",col="red", xlab="", ylab="")
legend(x=-0.08, y=60, legend = "calibrated Thallatosuchia",col = "red",lty=1, bty="n", y.intersp = 1.1)
axis(side=1,cex.axis=1.0, padj=1,at=seq(from=-(max(nodeHeights(tree.calibrated)+offset.calibrated)), to=0, by=(max(nodeHeights(tree.calibrated))+offset.calibrated)/23), labels=seq(230,0,by=-10))


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

###########################################################################
######################## Get the data #####################################
###########################################################################


### Load the tree from the working directory
Bayes.tree<-read.nexus("Data/Consensus_trees/MrBayes_calibrated_Thalattosuchia_201Ma_March2016.nex.con.tre")
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

### make an offset for the axis as R won't draw it from the tip to the root. The offset is a negative starting point for the axis equivalent to the
### rounding up we do at the root end of the axis i.e. if we round 227 Mya to 230 Mya then we need to offset by minus 3Ma of distance measured in
### branch lengths. To do this we divide the root height by the root age and multiply by -3.
offset<-3*(max(nodeHeights(tree.rename)/227))

### Read in the geological epoch data
epochs<-read.csv("Data/Metadata/epochs_and_colours.csv", stringsAsFactors = F)

### calculate the x co-ordinates for the calibrated tree polygons
epochs$calibrated.start<-(max(nodeHeights(tree.rename))-((max(nodeHeights(tree.rename)+offset)/23)*(epochs$Starting/10)))
epochs$calibrated.end<-(max(nodeHeights(tree.rename))-((max(nodeHeights(tree.rename)+offset)/23)*(epochs$Ending/10)))

### set up a vector of the epoch colours based on Commission for the Geological Map of the World guidelines
legend.cols<-rgb(red=epochs$Red, green=epochs$Green, blue=epochs$Blue, alpha=100, maxColorValue = 255)

## @knitr constrained_201_treeplot

#########################################################
### Plot the main combined tree with support values
#pdf(file=paste(out,file.names[[x]],"_",title,".pdf",sep=""), 30, 30)
plot(tree.rename,
     show.node.label=FALSE,
     cex=0.5,
     x.lim=0.5,
     label.offset=0.001)
nodelabels(tree.rename$node.label,adj=c(1,1),frame="none",
           col=ifelse(tree.rename$node.label>0.9,"red",
                      ifelse(tree.rename$node.label>=0.75 & tree.rename$node.label<0.9,"blue","#0000ff00")),cex=0.5)

## put on a the correct axis
axis(side=1,cex.axis=0.5,padj=1,at=seq(-offset,max(nodeHeights(tree.rename)),by=(max(nodeHeights(tree.rename))+offset)/23), labels=seq(230,0,by=-10))

### epoch plot order for trees
epoch.order<-seq(from=1, to=13, by=1)

### plot the epochs
for(i in epoch.order){
  polygon(x=c(epochs$calibrated.start[i],epochs$calibrated.start[i],epochs$calibrated.end[i],epochs$calibrated.end[i]),
          y = c(0,120,120,0), border=NA, col =legend.cols[i])
}
legend("topright", title="Epoch", inset=0.005, legend = epochs$Stage,
       fill =legend.cols,
       cex=0.5,
       bg = "white")


#dev.off()

# @knitr calibratedlttplot

### Lineages through time.
#pdf(file="Diagrams/calibrated_LTT_plot,pdf")
ltt.plot(tree.rename, xaxt="n", xlab="Time (Ma)", ylab="Extant lineages")
axis(side=1,cex.axis=1.0, padj=1,at=seq(from=-(max(nodeHeights(tree.rename)+offset)), to=0, by=(max(nodeHeights(tree.rename))+offset)/23), labels=seq(230,0,by=-10))
#dev.off()

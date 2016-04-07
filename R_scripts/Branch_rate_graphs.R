###########################################################################
############# Script for plotting rates of morphological change over time #################
###########################################################################

## @knitr constrained_201_treeprocess

### Clear the workspace
rm(list=ls())

# Load the Claddis package into R:
library(Claddis)
# Load the paleotree library into R (we will use this for time-scaling our tree later):
library(paleotree)
# Load the strap library into R (we will use this for plotting a time-scaled tree):
library(strap)

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

# You can download this file and read it in to R from your own hard drive, but for now we will read it in direct from the web as this means a single address (should!) work for everyone:
morpho.data <- ReadMorphNexus("Data/Metadata/Thalattosuchia_March2016_nexusmorpho.nex.txt")

# In the above line we have stored (<-) the data in a variable (nexus.data).
# This data is stored in an R format known as a "list", which is (one) way to store different types of data together.
# The components of this list (in this case) have names that can be seen by typing:
names(morpho.data)

# We can view individual components of this list by using the dollar symbol and a name from above, e.g.:
morpho.data$matrix

# Now we have data can try using some of the other functions.
# One function that we can run on just a morphological matrix (with no additional data) is Safe Taxonomic Reduction (Wilkinson 1995; Systematic Biology).
# This is a simple way to remove taxa we know (under parsimony) can only fall out in particular place(s) in the tree.
# We can thus "safely" remove them prior to inference and save ourselves some computation time.
# Here we will run it and store the results in a new variable (safe.data):
safe.data <- SafeTaxonomicReduction(morpho.data)

# This is also a list. For now we can just look at the first part of it:
safe.data$str.list

# You should see a matrix with three columns (Junior, Senior, and Rule).
# Here "Junior" is a taxon that can be safely removed, "Senior" is the taxon (or taxa) it will be expected to fall out next to and "Rule" is the rule (see Willinson 1995) under which it can be safely removed.
# In this case there is only one taxon ("Dry_Island_bonebed_specimens") that can be removed.
# Note that if you run this on your own data you may not be able to remove any taxa.
# If this occurs instead of a list you will simply get a string (text) telling you "No taxa can be safely removed".

# Another function we can use on just a matrix is MorphDistMatrix, which converts a cladistic matrix into a distance matrix using the various diferent metrics I mentioned in my talk: http://www.slideshare.net/graemelloyd/new-methodologies-for-the-use-of-cladistictype-matrices-to-measure-morphological-disparity-and-evolutionary-rate
# Again we will use a new variable (dist.data) to store the output:
dist.data <- MorphDistMatrix(morpho.data)

# This is (again!) a list. We can see the names of each part again using names():
names(morpho.data)

# These are the four metrics I discussed in my talk (max = MOD), plus an additional matrix giving the number of characters that are scored in BOTH taxa for each pairwise comparison.
# Lets check this last one quickly to see if there are any zeroes (i.e., incalculable distances):
any(dist.data$comp.char.matrix == 0)

# As this is true and most analyses don't appreciate incalculable distances, we can take a quick look at it to see if we can spot where the problem is (this time using the MOD metric):
dist.data$max.dist.matrix

# We can't use the "NAs" so we can use another function (TrimMorphDistMatrix) and another new variable (trimmed.max.data) that removes as few taxa as possible while dropping "NAs":
trimmed.max.data <-TrimMorphDistMatrix(dist.data$max.dist.matrix)

# This is (again!) a list.
# We can see what taxa have been removed by typing:
trimmed.max.data$removed.taxa

# We now also have a distance matrix without any gaps (NAs):
any(is.na(trimmed.max.data$dist.matrix))

# This should be FALSE, which means we can hand it to cmdscale and not get an error:
cmdscale(trimmed.max.data$dist.matrix)

# This should give a two-column matrix, but this is not what we really want as forcing all the variance on to two axes will confer a cost in terms of fidelity to the true distance.
# We can maximise our axes by upping the value "k" (an option in the function) to N - 1 (the maximum number of axes for N objects, i.e., N taxa).
# In addition we want to use another option in the function (add) which gets around the negative eigenvalue problem that can cause downstream problems (e.g., a scree plot with negative values).
# We can specify these options fairly easily and store our answer in a new variable (pco.data) and this time we will just use part of the output ($points) which are the values for our taxa on every ordination axis:
pco.data <- cmdscale(trimmed.max.data$dist.matrix, k=nrow(trimmed.max.data$dist.matrix) - 1, add=T)$points

# Before we plot this lets get the data we need to make a scree plot:
scree.data <- apply(pco.data, 2, var) / sum(apply(pco.data, 2, var)) * 100

# We can make a simple plot of this:
plot(scree.data, type="l", xlab="Ordination axis", ylab="Percentage variance")


# Before we start plotting a useful thing to do is define our plotting axes first so we can edit these later to easily plot different axes:
PCOx <- 1
PCOy <- 2

# We can use this line to plot our data:
plot(pco.data[, PCOx], pco.data[, PCOy], xlab=paste("PCO ", PCOx, " (", round(scree.data[PCOx], 2), "% variance)", sep=""), ylab=paste("PCO ", PCOy, " (", round(scree.data[PCOy], 2), "% variance)", sep=""), pch=19)

name.colours<-read.csv("Data/Metadata/Thallatosuchia_colours.csv", header=TRUE, stringsAsFactors = FALSE)

# This can be a bit tricky to interpret, not least of all because we do not know which point is which taxon.
# We can fix this by adding the taxon names:
text(pco.data[, PCOx], pco.data[, PCOy], rownames(pco.data), col= name.colours$colour[match(rownames(pco.data),name.colours$Species)])

# You can plot your tree in a much prettier way using Mark Bells awesome geoscalePhylo function in the strap package:
geoscalePhylo(ladderize(My.tree), cex.age=0.6, cex.ts=0.8, cex.tip=1)

### remove polytomies - we might want to do this with MCC trees in the end.
My.tree <- multi2di(My.tree)

### give zero length branches an arbitrarily short length
My.tree$edge.length<-ifelse(My.tree$edge.length==0,0.001,My.tree$edge.length)

# Now we have a time-scaled tree we can calculate rates.
# As I (briefly) mentioned this is hard to do for time series so please disregard that output from the function (hopefully I will get a good approach working in future though).
# At the moment the function still wants some time bins so in the below we will simply set six equally spaced time bins from the root to the youngest tip.
# Another option in the function is to set the alpha value for the significance tests.
# Note that the function already accounts for multiple comparisons using the Benjamini and Hochberg 1995 JRSSB False Discovery Rate (FDR).
# We can run the rate tests using DiscreteCharacterRate and store the results in a new variable (rate.data):
rate.data <- DiscreteCharacterRate(My.tree, morpho.data, seq(My.tree$root.time, My.tree$root.time - max(diag(vcv(My.tree))), length.out=6), alpha=0.01)

# Again this is a list and the output is quite verbose, e.g....:
rate.data$branch.results

# However, a better way to look at the data is visually by colouring the branches.
# We can start by creating a vector of colours for our branches:
edge.color <- rep("black", nrow(My.tree$edge))
edge.color[which(rate.data$branch.results[, "ml.signif.hi"] == 1)] <- "red"
edge.color[which(rate.data$branch.results[, "ml.signif.lo"] == 1)] <- "blue"

pdf("Diagrams/relative_rates.pdf",30,20)
# We can now plot our tree with branches coloured by rate (black = non-significant rates, red = significantly high rates, blue = significantly low rates):
geoscalePhylo(ladderize(My.tree), cex.age=0.6, cex.ts=0.8, cex.tip=1, edge.color=edge.color[match(ladderize(My.tree)$edge[, 2], My.tree$edge[,2])])
dev.off()
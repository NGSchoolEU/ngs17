source("scripts/plotTree.R")
library(phytools)
library(phangorn)

tru <- read.tree("phylogeny/RAxML_bipartitions.NC_000962.b1")
trr <- midpoint(tru)

write.tree(trr, file="temp.tre")
trr2 <- read.tree("temp.tre")

plotTree(trr2, infoFile="scripts/tipinfo.csv", colourNodesBy ="case", tip.colour.cex=1, lwd=1, tipColours = c("blue","red", "green", "black","purple","purple2","skyblue2","grey"), tip.labels=T, offset=0.03, tipLabelSize=1, legend=F, infoCols=NA, heatmapData="scripts/resist.csv",heatmap.colours=c("black","lightgrey"), treeWidth=10,dataWidth=2, colLabelCex=1.5, outputPNG="tree.png",w=1000,h=1500,closeDev=F)


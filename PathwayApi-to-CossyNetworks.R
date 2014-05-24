#####################################################
## this script creates MISs from Cossy Network. #####
## Cossy networks has 3 types of files. 
## One file contains the names of pathways (subnetworks), FileName: pathway-***.txt
## Node files (FileName: node-[PATHWAY].txt) contains a map from node id to genes.
## Edge files (FileName: edge-[PATHWAY].txt) contains the edges
## all generated MISs are saved in a gmt file.
#####################################################

## settings
pathwayApiDataDir <- 'data/pathwaycsv'
network <- 'pathwayapi'
outputNetworkDir <- 'results'

pathwayApiRelationsFile <- paste0(pathwayApiDataDir, '/pathway_rs.csv')
geneMappingFile <- paste0(pathwayApiDataDir, '/gene_mapping.csv')


## read input file
geneRelations <- read.csv(file=pathwayApiRelationsFile, header=T, sep=";", quote="\"",  stringsAsFactors=F, comment.char="")
geneMappings <- read.csv(file=geneMappingFile, header=T, sep=";", quote="\"",  stringsAsFactors=F, comment.char="", row.names=1)

## check if any gene name contains a white space
if(length(grep("\\s", geneMappings[,'gene_name']))>0)
  stop("gene name cannot have any white space")



## Take only gene interactions
geneRelations <- geneRelations[,2:3,drop=F]

## convert gene_id to gene_names in gene-relations
getGeneName <- function(geneid){
  return(geneMappings[geneid,])
}
geneRelations[,1] = getGeneName(geneRelations[,1])
geneRelations[,2] = getGeneName(geneRelations[,2])

## remove any relation if it does not have any name
validRelations <- !is.na(geneRelations[,1]) & !is.na(geneRelations[,2])
geneRelations <- geneRelations[validRelations,,drop=F]

## Each unique gene will be a node
uniqueGenes <- unique(c(geneRelations[,1], geneRelations[,2]))

## give a nodeID to each unique gene, NodeID must start from 1
nodeIDs <- setNames(1:length(uniqueGenes), uniqueGenes)

## convert gene_names to node_id in gene-relations
getNodeID <- function(geneName){
  return(nodeIDs[geneName])
}
geneRelations[,1] = getNodeID(geneRelations[,1])
geneRelations[,2] = getNodeID(geneRelations[,2])



## write output files - 1) network file
netFile <- paste0(outputNetworkDir, "/", network, ".txt")
cat(paste0(network, "\n"), file=netFile)

## write output files - 2) node files
nodeFile <- paste0(outputNetworkDir, "/node-", network, ".txt")
nodeIDsDf <- data.frame(id=nodeIDs, name=names(nodeIDs))
write.table(nodeIDsDf, file=nodeFile, quote=F, sep="\t", row.names=F, col.names=F)

## write output files - 3) edge files
edgeFile <- paste0(outputNetworkDir, "/edge-", network, ".txt")
write.table(geneRelations, file=edgeFile, quote=F, sep="\t", row.names=F, col.names=F)

## Show notifications
cat(paste0("Generated 3 files:\n1. ", netFile, "\n2. ", nodeFile, "\n3. ", edgeFile))


#############################################################################################
## this script creates CossyNetwork (node and edge files) from original kegg files #####
## Cossy networks has 3 types of files.
## One file contains the names of pathways (subnetworks), FilePAth: data/[NETWORK]/[NETWORK].txt
## Node files (FileName: node-[PATHWAY].txt) contains a map from node id to genes.
## Edge files (FileName: edge-[PATHWAY].txt) contains the edges
## The generated outputs are saved in 'results' folder
############################################################################################

library(biomaRt)
library(plyr)
library(graph)
library(KEGGgraph)

## input files
kgmlDir <- 'data/raw/kegg-2014-07-16/'

## output files
outPathwaysFile <- 'results/kegg.txt'
outNodeFile <- 'results/node-kegg.txt'
outEdgeFile <- 'results/edge-kegg.txt'

## settings
invalid_gene_id <- "2147483647"

#### function to get gene symbols from gene ids ########
geneSymbolsFromGeneIds <- function(gids, martobj){
  
  if(length(unique(gids))!=length(gids))
    stop("Every gene id must be unique.")
  
  gid2sym = getBM(attributes=c('entrezgene','hgnc_symbol'), filters='entrezgene', values=gids, mart=martobj)
  
  symbols <- rep(NA, length(gids))
  names(symbols) <- gids
  # multiple entries may exist. Use the longest value to avoid blank entry.
  symbolsFromBiomart <- tapply(gid2sym$hgnc_symbol, gid2sym$entrezgene, max)
  symbols[names(symbolsFromBiomart)] = symbolsFromBiomart
  symbols[symbols==""] <- NA
  
  return(symbols)
}

getInteractionsFromKgml <- function(kgmlFile){
  ## create an empty dataframe to return in any error case
  emptyDF <- data.frame(from="dummy", to="dummy", stringsAsFactors = F)[-1,,drop=F]
  
  ## create a graph
  kgraph <- try(parseKGML2Graph(kgmlFile, expandGenes=TRUE), silent = T)
  if(class(kgraph) == "try-error")
    return(emptyDF)
  
  ## get edges
  kedges <- edges(kgraph)
  
  ## create a data frame from edges
  len <- length(kedges)
  from <- c()
  if(len > 1){
    for(fr in names(kedges))
      from <- c(from, rep(fr, length(kedges[[fr]])))
  }
  to <- unlist(kedges, use.names = F)
  
  if(length(from)==0)
    return(emptyDF)
  
  edgeDF <- data.frame(from=from, to=to, stringsAsFactors = F)
  return(edgeDF)
}

## read files
kgmlFiles <- list.files(path = kgmlDir, pattern = "*\\.kgml$")
kgmlFiles <- paste0(rep(kgmlDir, length(kgmlFiles)), kgmlFiles)
allEdgeFrames <- lapply(kgmlFiles, getInteractionsFromKgml)
edgeDF <- ldply(allEdgeFrames, data.frame)


###### create pathways file. ########
cat("kegg\n", file = outPathwaysFile, sep = "")


###### create node file  #########

## get all unique proteins and corresponding gene ids and symbols
keggGenes <- unique(c(edgeDF$from, edgeDF$to))
geneids <- sapply(keggGenes, function(kid) strsplit(kid,"hsa:")[[1]][2])
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
geneidToSymbols <- geneSymbolsFromGeneIds(geneids, ensembl)


## write in node file
nodeMap <- data.frame(id=1:length(geneids), geneid=geneids, genesymbol=geneidToSymbols, row.names = keggGenes)
write.table(nodeMap[,c("id","genesymbol")], file = outNodeFile, quote = F, sep = "\t", na = invalid_gene_id, row.names = F, col.names = F)

####### create edge file #########
fromNodeIds <- sapply(edgeDF$from, function(fr) nodeMap[fr,"id"])
toNodeIds <- sapply(edgeDF$to, function(fr) nodeMap[fr,"id"])
netEdgeDF <- data.frame(from=as.integer(fromNodeIds), to=as.integer(toNodeIds))

## make undirectional and remove multiple edges
from <- apply(netEdgeDF, 1, min)
to <- apply(netEdgeDF, 1, max)
netEdgeDF <- data.frame(from=from, to=to)
netEdgeDF <- unique(netEdgeDF)

# writing edges
write.table(netEdgeDF, file = outEdgeFile, quote = F, sep = "\t", row.names = F, col.names = F)

####### show notification ################
cat("CossyNetwork saved in results folder. 1)", outPathwaysFile, " 2)", outNodeFile, " 3)", outEdgeFile)

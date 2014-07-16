#############################################################################################
## this script creates CossyNetwork (node and edge files) from original string db files #####
## Cossy networks has 3 types of files.
## One file contains the names of pathways (subnetworks), FilePAth: data/[NETWORK]/[NETWORK].txt
## Node files (FileName: node-[PATHWAY].txt) contains a map from node id to genes.
## Edge files (FileName: edge-[PATHWAY].txt) contains the edges
## The generated outputs are saved in 'results' folder
############################################################################################

library(biomaRt)
library(plyr)

## input files
stringInteractionFile <- 'data/raw/string-v9.1/9606.protein.links.v9.1.txt'
geneMapFile <- 'data/raw/string-v9.1/entrez_gene_id.vs.string.v9.05.28122012.txt'

## output files
outPathwaysFile <- 'results/string.txt'
outNodeFile <- 'results/node-string.txt'
outEdgeFile <- 'results/edge-string.txt'

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

## read files
interactions <- read.table(file = stringInteractionFile, header = T, sep = "", quote = "", colClasses = c('character', 'character', 'numeric'), col.names = c("protein1","protein2","score"), comment.char = "")
geneMaps <- read.table(file = geneMapFile, header = T, sep = "", quote = "", colClasses = c('character', 'character'), col.names = c("geneid","proteinid"), comment.char = "")

## take only high confidence interactions (score > 0.7)
interactions <- interactions[interactions$score>700, ]

###### create pathways file. ########
cat("string\n", file = outPathwaysFile, sep = "")


###### create node file  #########

## get all unique proteins and corresponding gene ids and symbols
proteins <- unique(c(interactions$protein1, interactions$protein2))
proteinidToGeneidMap <- tapply(geneMaps$geneid, INDEX = geneMaps$proteinid, function(pid) pid)

## Some proteins (locus id) do not have gene mapping (strange, but true! For e.g. 9606.ENSP00000382366 in v9.1)
## so, add dummy invalid geneid to proteinidToGeneidMap to maintain the network information.
unmappedProteins <- setdiff(proteins, geneMaps$proteinid)
if(length(unmappedProteins) > 0){
  for(up in 1:length(unmappedProteins)){
    unmappedProtien <- unmappedProteins[up]
    proteinidToGeneidMap[[unmappedProtien]] <- paste0("UnmappedGene", up)
  }
}

geneids <- unique(unlist(proteinidToGeneidMap))
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
geneidToSymbols <- geneSymbolsFromGeneIds(geneids, ensembl)


## write in node file
nodeMap <- data.frame(id=1:length(geneids), geneid=geneids, genesymbol=geneidToSymbols, row.names = geneids)
write.table(nodeMap[,c("id","genesymbol")], file = outNodeFile, quote = F, sep = "\t", na = invalid_gene_id, row.names = F, col.names = F)

####### create edge file #########

emptyEdgeDataFrame <- data.frame(from="dummy", to="dummy", stringsAsFactors = F)[-1,,drop=F]

constructEdges <- function(proteinInteraction){
  p1 <- proteinInteraction[1]
  p2 <- proteinInteraction[2]
  
  # proteinid to geneids
  geneids1 <- proteinidToGeneidMap[[as.character(p1)]]
  geneids2 <- proteinidToGeneidMap[[as.character(p2)]]
  
  # geneid to nodeid
  nodeids1 <- nodeMap[geneids1, "id"]
  nodeids2 <- nodeMap[geneids2, "id"]
  
  if(length(nodeids1) == 0 || length(nodeids2) == 0)
    print(proteinInteraction)
  
  lenN1 <- length(nodeids1)
  lenN2 <- length(nodeids2)
  
  from <- c()
  for(n1 in nodeids1){
    from <- c(from, rep(n1, lenN2))
  }
  to <- rep(nodeids2, lenN1)
  edges <- data.frame(from=from, to=to, stringsAsFactors = F)
  
  return(edges)
}

## construct edges
allEdges <- apply(interactions, 1, constructEdges)

## make an edge dataframe
edgeDF <- ldply(allEdges, data.frame)
edgeDF <- edgeDF[,c("from","to")]           # take only from and to columns, NOT id colums
edgeDF$from <- as.integer(edgeDF$from)      # make integer
edgeDF$to <- as.integer(edgeDF$to)          # make integer

## make undirectional and remove multiple edges
from <- apply(edgeDF, 1, min)
to <- apply(edgeDF, 1, max)
edgeDF <- data.frame(from=from, to=to)
edgeDF <- unique(edgeDF)

#print("writing edges")
write.table(edgeDF, file = outEdgeFile, quote = F, sep = "\t", row.names = F, col.names = F)


############## V-0
# writeInteractionToEdgeFile <- function(proteinInteraction, outFile){
#   p1 <- proteinInteraction[1]
#   p2 <- proteinInteraction[2]
#   
#   # proteinid to geneids
#   geneids1 <- proteinidToGeneidMap[[as.character(p1)]]
#   geneids2 <- proteinidToGeneidMap[[as.character(p2)]]
#   
#   # geneid to nodeid
#   nodeids1 <- nodeMap[geneids1, "id"]
#   nodeids2 <- nodeMap[geneids2, "id"]
#   
#   if(length(nodeids1) == 0 || length(nodeids2) == 0)
#     print(proteinInteraction)
#   
#   lapply(nodeids1, function(n1){
#     lapply(nodeids2, function(n2){
#       cat( paste0(n1, '\t', n2, '\n'), append = T, file = outFile)
#       return(NA)
#     })
#     return(NA)
#   })
#   
# }
# 
# 
# outEdgeConn <- file(outEdgeFile, "w")
# tmp <- apply(interactions, 1, writeInteractionToEdgeFile, outEdgeConn)
# close(outEdgeConn)

##################################

####### show notification ################
cat("CossyNetwork saved in results folder. 1)", outPathwaysFile, " 2)", outNodeFile, " 3)", outEdgeFile)

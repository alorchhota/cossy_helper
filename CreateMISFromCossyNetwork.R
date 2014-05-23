#####################################################
## this script creates MISs from Cossy Network. #####
## Cossy networks has 3 types of files. 
## One file contains the names of pathways (subnetworks), FileName: pathway-***.txt
## Node files (FileName: node-[PATHWAY].txt) contains a map from node id to genes.
## Edge files (FileName: edge-[PATHWAY].txt) contains the edges
## all generated MISs are saved in a gmt file.
#####################################################

library(igraph)

## settings
minNodeInGis <- 5
maxNodeInGis <- 25
directed <- F
approximate <- T
rmvMultipleEdges <- T      # it is required for generating mis from kegg network, not for string.

invalid_gene_id <- "2147483647"

## input-output files
pathwaysDataDir <- "data/networks/kegg"
pathwaysFileName <- "pathways-kegg.txt"
gmtOutFile <- paste("results/", pathwaysFileName, "-mis-", minNodeInGis, "-", maxNodeInGis, ".gmt", sep="")




bigGisSizes <- c()

bigGisLines <- c()
smallGisLines <- c()
okGisLines <- c()

extractGisFromKgml <- function (kgmlNodeFile, kgmlEdgeFile, rmvMultipleEdges=F){
  edgeData <- read.table(kgmlEdgeFile, col.names=c("from","to"))
  edgeData <- edgeData+1
  
  ## TODO: what to return
  if(nrow(edgeData) <=0 )
    return (-1) 
  
  ## remove multiple edges for undirected graph
  if(rmvMultipleEdges)
    edgeData <- removeMultipleEdges(edgeData, directed)
  
  
  
  nodeData <- read.table(kgmlNodeFile, sep="\t", col.names=c("id","genes"), colClasses=c("numeric","character"))
  nodeData$id <- nodeData$id+1
  
  kgmlGraph <- buildIGraph(nodeData, edgeData, directed)
  components <- clusters(graph=kgmlGraph)
  
  # for each component, extract GIS
  for(comp in 1:components$no){
    
#     if(comp==3)
#       debug(extractGisFromConnectedGraph)
    
    componentNodes <- which(components$membership==comp)
    componentGraph <- induced.subgraph(graph=kgmlGraph, vids=componentNodes)
    
    componentNodeData <- nodeData[componentNodes, ,drop=F]
    componentNodeData$id <- 1:nrow(componentNodeData)
    componentEdgeData <- get.edgelist(graph=componentGraph, names=F)
    
    gisLines <- extractGisFromConnectedGraph(componentNodeData, componentEdgeData)
    saveGisLines(gisLines)
    
    ##count low size network
    if(nrow(componentNodeData) < minNodeInGis){
      smallGisLines <<- c(smallGisLines, gisLines)
    }
    else if(nrow(componentNodeData) <= maxNodeInGis){
      okGisLines <<- c(okGisLines, gisLines)
    }
  }
  
}
#debug(extractGisFromKgml)

removeMultipleEdges <- function(edgeData, directed=FALSE){
  N <- nrow(edgeData)
  if(N<=1)
    return(edgeData)
  
  toRemove <- rep(FALSE, times=N)
  for(i in 1:(N-1)){
    
    if(toRemove[i]){
      next
    }
    
    ni1 <- edgeData[i,1]
    ni2 <- edgeData[i,2]
    
    for(j in (i+1):N){
      nj1 <- edgeData[j,1]
      nj2 <- edgeData[j,2]
      
      if(ni1==nj1 && ni2==nj2)
        toRemove[j] = TRUE
      
      if(directed==FALSE && ni1==nj2 && ni2==nj1)
        toRemove[j] = TRUE
    }
  }
  
  return(edgeData[!toRemove,])
}

extractGisFromConnectedGraph <- function(nodeData, edgeData){
  if(nrow(nodeData)==0){
    return("")
  }
  
  if(nrow(nodeData)<minNodeInGis){
    # taking gis from network having less than min nodes
    return(genesInGraph(nodeData))
  }
  
  if(nrow(nodeData)<=maxNodeInGis){
    # the network is in range, take the full network
    return(genesInGraph(nodeData))
  }
  
  kgmlGraph <- buildIGraph(nodeData, edgeData, directed)
  noOfNodes <- vcount(kgmlGraph)
  
  if(approximate)
    bet <- fastgreedy.community(kgmlGraph)
  else
    bet <- edge.betweenness.community(kgmlGraph)
  
  
  graphList <- apply(nodeData, 1, function(node){
    g <- graph.empty(n=0, directed=F)
    g <- set.graph.attribute(graph=g, name="genes", value=node["genes"])
    g <- set.graph.attribute(graph=g, name="count", value=1)        # 1 node, count=1
    g <- set.graph.attribute(graph=g, name="nodes", value=as.integer(node["id"]))
    #g <- set.graph.attribute(graph=g, name="edge", value="")
    return(g)
  })
  
  # merge
  graphListLength <- length(graphList)
  mergeWithGraphList <- function(mergeIndexes){
    g1 <- graphList[[mergeIndexes[1]]]
    g2 <- graphList[[mergeIndexes[2]]]
    g <- mergeGraphs(g1,g2)
    graphListLength <<- graphListLength+1
    graphList[[graphListLength]] <<- g
    return(NA)
  }
  apply(bet$merges, 1, mergeWithGraphList)
  
  # find community
  isTaken <- rep(FALSE, graphListLength)
  takenGraphs <- list()
  discardedGraphs <- list()
  takenGraphsCount <- 0
  discardedGraphsCount <- 0
  
  findCommunityGraphs <- function(mergeRowIndex){
    graphIndex <- noOfNodes+mergeRowIndex
    index1 <- bet$merges[mergeRowIndex,1]
    index2 <- bet$merges[mergeRowIndex,2]
    
    g <- graphList[[graphIndex]]
    
    if(isTaken[graphIndex]){
      # the graph is taken means, its subgraphs are also taken
      isTaken[index1] <<- TRUE
      isTaken[index2] <<- TRUE
    }
    else {
      g.count <- as.integer(get.graph.attribute(graph=g, name="count"))
      if(g.count<=maxNodeInGis){
        if(g.count>=minNodeInGis){
          takenGraphs[[takenGraphsCount+1]] <<- g
          takenGraphsCount <<- takenGraphsCount+1
        }
        else{
          discardedGraphs[[discardedGraphsCount+1]] <<- g
          discardedGraphsCount <<- discardedGraphsCount+1
        }
        
        isTaken[index1] <<- TRUE
        isTaken[index2] <<- TRUE
      }
    }
    
    isTaken[graphIndex] <<- TRUE
    
    return(NA)
  }
  
  lapply(nrow(bet$merges):1, findCommunityGraphs)
  
  
  # take discarded single nodes
  takeDiscardedSingleNodes <- function(nodeIndex){
    if(!isTaken[nodeIndex]){
      discardedGraphs[[discardedGraphsCount+1]] <<- graphList[[nodeIndex]]
      discardedGraphsCount <<- discardedGraphsCount+1
    }
    return(NA)
  }
  lapply(1:noOfNodes, takeDiscardedSingleNodes)
  
  
#   if(length(takenGraphs)==0 && length(discardedGraphs)==1){
#     # The graph has less than minNodeInGis
#     # prepare a gis and return it
#     gisGenes <- genesInAnnotatedGraph(discardedGraphs[[1]])
#     return(gisGenes)
#   }
#   
#   if(length(takenGraphs)==1 && length(discardedGraphs)==0){
#     # all genes are taken in one graph i.e. give graph size is within the range
#     gisGenes <- genesInAnnotatedGraph(takenGraphs[[1]])
#     return(gisGenes)
#   }
  
  
  # find the nodes within 'among' that have distance d from n
  assignNode <- function(n, among, d){
    #candidates <- nodes in among that have d distance
    candidates <- getDHopNodes(n, among,d)
    
    # return if no d-distant nodes found
    if(length(candidates) == 0){
      if(d==1){
        # no connection with any gis => no membership
        return(-1)
      }
      else{
        # connected with few gis, but could not decide from less distance
        # each GIS is equal candidate
        # should return a random gis (the gis of any random node)
        candidates <- sample(among,1)
      }
    }
    
    gisFreq <- countGisMembershipFrequency(candidates)
    maxFreqGis <- which(gisFreq==max(gisFreq))
    
    # if(only 1 gis have the highest fequency)
    #   return that gis
    if(length(maxFreqGis)==1){
      return(maxFreqGis)
    }
    else{
      # if (more than 1 gis have the highest frequency)
      #   nodes <- vector of nodes of the giss with highest feq
      #   assignNode(n, nodes, d+1)
      
      candidateGisGraphs <- lapply(maxFreqGis, function(gisIndex){return(takenGraphs[[gisIndex]])})
      newAmong <- getAllNodes(candidateGisGraphs)
      newAssignment <- assignNode(n, newAmong, d+1)
      return(newAssignment)
      
    }
    
    return(-1)  
  }
  #debug(assignNode)
  
  getDHopNodes <- function(from, among, d){
    dNeihbors <- neighborhood(graph=kgmlGraph, nodes=from, order=d)[[1]]
    dLess1Neihbors <- neighborhood(graph=kgmlGraph, nodes=from, order=d-1)[[1]]
    onlyDNeighbors <- setdiff(dNeihbors, dLess1Neihbors)
    
    targetDNeighbors <- intersect(onlyDNeighbors, among)
    
    # remove unused variables
    rm(dNeihbors, dLess1Neihbors, onlyDNeighbors)
    
    return(targetDNeighbors)
  }
  
  countGisMembershipFrequency <- function(nodes){
    #       gisFrequency <- 0 for all gis
    #       for each c in candidates
    #         gis <- the gis c is in
    #         increase gisFrequency for the gis
    
    freq <- rep(0, length(takenGraphs))
    lapply(nodes, function(n){
      gi <- findConnectingGraphIndex(graphs=takenGraphs, n)
      freq[gi] <<- freq[gi] + 1
    })
    
    return(freq)
  }
  
#   addNodesInGisGraphs <- function(discardedNodes, memberships){
#     discardedNodeIndex <- 0
#     addedNodes <- c()
#     lapply(discardedNodes, function(dn){
#       discardedNodeIndex <<- discardedNodeIndex + 1
#       membershipGraphIndex <- memberships[discardedNodeIndex]
#       nodeIndex <- discardedNodes[discardedNodeIndex]
#       if(membershipGraphIndex > 0){
#         takenGraphs[[membershipGraphIndex]] <<- mergeGraphs(takenGraphs[[membershipGraphIndex]], graphList[[nodeIndex]])
#         addedNodes[length(addedNodes)+1] <<- nodeIndex
#       }
#     })
#     
#     return(addedNodes)
#   }
  
  # if any node or community is discarded, assign them to some taken graph.
  if(length(discardedGraphs) > 0){
    takenNodes <- getAllNodes(takenGraphs)
    discardedNodes <- getAllNodes(discardedGraphs)
    
    while(length(discardedNodes) > 0){
      print(paste(format(Sys.time(), "%X"), "#discarded nodes:",length(discardedNodes)))
      memberships <- sapply(discardedNodes, assignNode, takenNodes, 1)
      
#       testDebug <- function(){
#         tmp <- 0
#       }
#       debug(testDebug)
#       testDebug()
      
      if(sum(memberships>0) == 0){
        # no new memberships, so break merging
        break
      }
      
      #newlyAddedNodes <- addNodesInGisGraphs(discardedNodes, memberships)
      addNodesInGisGraphs <- function(){
        discardedNodeIndex <- 0
        addedNodes <- c()
        lapply(discardedNodes, function(dn){
          discardedNodeIndex <<- discardedNodeIndex + 1
          membershipGraphIndex <- memberships[discardedNodeIndex]
          nodeIndex <- discardedNodes[discardedNodeIndex]
          if(membershipGraphIndex > 0){
            takenGraphs[[membershipGraphIndex]] <<- mergeGraphs(takenGraphs[[membershipGraphIndex]], graphList[[nodeIndex]])
            addedNodes[length(addedNodes)+1] <<- nodeIndex
          }
        })
        
        return(addedNodes)
      }
      newlyAddedNodes <- addNodesInGisGraphs()
      
      
      # update discarded and taken node vector
      #discardedNodes <- discardedNodes[-which(discardedNodes %in% newlyAddedNodes)]
      discardedNodes <- setdiff(discardedNodes, newlyAddedNodes)
      takenNodes <- c(takenNodes, newlyAddedNodes)
      
    }
    print(paste(format(Sys.time(), "%X"), "#discarded nodes:",length(discardedNodes)))
    
  }
  
  if(takenGraphsCount==1){
    # could not divide more, so return it. Or divide in special way
    #print("Could not resplit")
    bigGisSizes[length(bigGisSizes)+1] <<- get.graph.attribute(graph=takenGraphs[[1]], name="count")
    gisGenes <- genesInAnnotatedGraph(takenGraphs[[1]])
    
    ##count
    bigGisLines <<- c(bigGisLines, gisGenes)
    
    return(gisGenes)
  }
  else{
    # gisGenes
    # for every taken graph
    #   if(graph size <= maxThreshold)
    #     gisGenes.add(gisFromGraphs)
    #   else  
    #     create new node and edge data and recursively call this method
    #     gisGenes.add(result from above method)
    
    gisGenes <- c()
    lapply(takenGraphs, function(tg){
      tg.count <- as.integer(get.graph.attribute(graph=tg, name="count"))
      if(tg.count<=maxNodeInGis){
        gline <- genesInAnnotatedGraph(tg)
        gisGenes <<- c(gisGenes, gline)
        
        ## count
        okGisLines[length(okGisLines)+1] <<- gline
      }
      else{
        
        print(paste(format(Sys.time(), "%X"), "re-split:", tg.count))
        
        nodeIndexes <- getIntegerNodeVector(tg)
        nodeIndexes <- sort(x=nodeIndexes, decreasing=F)
        
        tgGraph <- induced.subgraph(graph=kgmlGraph, vids=nodeIndexes)
        tgNodeData <- nodeData[nodeIndexes, ,drop=F]
        tgNodeData$id <- 1:nrow(tgNodeData)
        tgEdgeData <- get.edgelist(graph=tgGraph, names=F)
        
        tgGisGenes <- extractGisFromConnectedGraph(tgNodeData, tgEdgeData)
        #print(tgGisGenes)
        gisGenes <<- c(gisGenes, tgGisGenes)
      }
    })
    
    return(gisGenes)
  }
}





buildIGraph <- function(nodeData, edgeData, directed=F){
  noOfNodes <- nrow(nodeData)
  edgeMatrix <- as.matrix(edgeData)
  edgeVector <- as.vector(t(edgeMatrix))
  gr <- graph(edges=edgeVector, n=noOfNodes, directed=directed)
  return(gr)
}

mergeGraphs <- function(g1, g2){
  
  g1.nodes <- get.graph.attribute(graph=g1, name="nodes")
  g1.genes <- get.graph.attribute(graph=g1, name="genes")
  g1.count <- as.integer(get.graph.attribute(graph=g1, name="count"))
  
  g2.nodes <- get.graph.attribute(graph=g2, name="nodes")
  g2.genes <- get.graph.attribute(graph=g2, name="genes")
  g2.count <- as.integer(get.graph.attribute(graph=g2, name="count"))
  
  g <- graph.empty()
  separator <- ifelse(nchar(g1.nodes) > 0 && nchar(g2.nodes)>0, ",", "")
  g <- set.graph.attribute(graph=g, name="nodes", value=paste(g1.nodes, separator, g2.nodes, sep=""))
  separator <- ifelse(nchar(g1.genes) > 0 && nchar(g2.genes)>0, ",", "")
  g <- set.graph.attribute(graph=g, name="genes", value=paste(g1.genes, separator, g2.genes, sep=""))
  g <- set.graph.attribute(graph=g, name="count", value=g1.count+g2.count)
  
  return(g)
}

genesInAnnotatedGraph <- function(g){
  g.genes <- get.graph.attribute(graph=g, name="genes")
  return(getValidSortedUniqueGenes(g.genes))
}

genesInGraph <- function(nodeData){
  genes <- ""
  lapply(nodeData$genes, function(g){
    if(genes=="")
      genes <<- g
    else
      genes <<- paste(genes,g,sep=",")
  })
  return(getValidSortedUniqueGenes(genes))
}

trim <- function(str){
  str <- gsub(pattern="^[ \t]*", replacement="", x=str)
  str <- gsub(pattern="*[ \t]*$", replacement="", x=str)
  str
}

getValidSortedUniqueGenes <- function(g.genes){
  g.genes <- strsplit(x=g.genes, split=",")[[1]]
  g.genes <- sort(unique(g.genes))
  g.genes <- setdiff(g.genes, invalid_gene_id)
  
  if(length(g.genes) == 0){
    return("")
  }
  
  line <- ""
  lapply(g.genes, function(gid){
    gid <- trim(gid)
    if(line=="") 
      line<<-gid
    else
      line <<- paste(line, gid, sep="\t")
  })
  
  return(line)
}

findConnectingGraphIndex <- function(graphs, n){
  index <- 0
  for(g in graphs){
    index <- index + 1
    g.nodes <- get.graph.attribute(graph=g, name="nodes")
    if(containNode(g.nodes, n)){
      return(index)
    }
  }
  
  return(-1)
}


containNode <- function(nodes, n){
  nodes <- as.character(nodes)
  splitted <- strsplit(nodes, split=',')[[1]]
  return(length(which(splitted==n)) > 0 )
}

getIntegerNodeVector <- function(g){
  nodes <- as.character(get.graph.attribute(graph=g, name="nodes"))
  splitted <- strsplit(nodes, split=',')[[1]]
  intnodes <- sapply(splitted, as.integer )
  return(intnodes)
}


getAllNodes <- function(graphList){
  nodes <- c()
  for(g in graphList){
    g.nodes <- getIntegerNodeVector(g)
    nodes <- c(nodes, g.nodes)
  }
  return(nodes)
}


allGis <- data.frame(id=c(), pathway=c(), genes=c(), stringsAsFactors=F)
saveGisLines <- function(gisGeneLines){
  
  lapply(gisGeneLines, function(line){
    if(line=="")
      return(NA)
    
    prevIndex <- which(allGis$genes == line)
    if(length(prevIndex)>0){
      if(!(allGis$pathway[prevIndex] == gmtPathway)){
        allGis$pathway[prevIndex] <<- paste(allGis$pathway[prevIndex], gmtPathway, sep="|")
      }
      
      #print(allGis)
    }
    else{
      newGis <- data.frame(id=gmtCounter, pathway=gmtPathway, genes=line, stringsAsFactors=F)
      allGis <<- rbind(allGis, newGis)
      gmtCounter <<- gmtCounter+1
    }
    
    return(NA)
  })
  
  return(NA)
}

printInGmtFile <- function(){
  
  apply(allGis, 1, function(gmtRow){
    line <- paste(gmtRow["id"], gmtRow["pathway"], gmtRow["genes"], sep="\t")
    write(file=gmtOut, x=line)  
  })
  
}


#debug(extractGisFromConnectedGraph)

### Extraction process ############

#read dataset
pathwaysFile <- paste0(pathwaysDataDir, "/", pathwaysFileName)
pathways <- read.table(pathwaysFile)

#### debug test for one pathway
#pathways <- pathways[229,,drop=F]
#pathways <- pathways[219,,drop=F]


## create gmt output file
gmtOut <- file(gmtOutFile, open="w+")
gmtCounter <- 0
gmtPathway <- "dummy"

#debug(saveGis)
apply(pathways, 1, function(pathwayName){
  kgmlNodeFile <- paste(pathwaysDataDir, "/node-", pathwayName, ".txt", sep="")
  kgmlEdgeFile <- paste(pathwaysDataDir, "/edge-", pathwayName, ".txt", sep="")
  gmtPathway <<- pathwayName
  
  #print(c(gmtPathway, nrow(allGis)) )
  print(gmtPathway)
  extractGisFromKgml(kgmlNodeFile, kgmlEdgeFile, rmvMultipleEdges=rmvMultipleEdges)
  
})

#saveGisLines(okGisLines)
nokgis <- nrow(allGis)
print(paste("ok gis:", nokgis))

#saveGisLines(bigGisLines)
nbiggis <- nrow(allGis) - nokgis
print(paste("big gis:", nbiggis))

#saveGisLines(smallGisLines)
nsmallgis <- nrow(allGis) - nokgis - nbiggis
print(paste("small gis:", nsmallgis))



printInGmtFile()
cat(gmtOutFile, "\n")
close(gmtOut)

cat("#GIS(Total):", nrow(allGis), "\n")
cat("#GIS(Big):", length(bigGisSizes), "\n")
cat("Big Sizes:", bigGisSizes, sep="\t")

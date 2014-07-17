#####################################################
## this script creates MISs from Cossy Network. #####
## Cossy networks has 3 types of files. 
## One file contains the names of pathways (subnetworks), FilePAth: data/[NETWORK]/[NETWORK].txt
## Node files (FileName: node-[PATHWAY].txt) contains a map from node id to genes.
## Edge files (FileName: edge-[PATHWAY].txt) contains the edges
## all generated MISs are saved in a gmt file.
#####################################################

library(igraph)
library(stringr)

source('cytoscape-network-json.R')

## settings
mis_init <- function(){
  allowedNetworks <<- c("kegg", "string", "pathwayapi")
  network <<- "kegg"
  minNodeInGis <<- 5
  maxNodeInGis <<- 15
  
  directed <<- FALSE     # molecular interaction networks are undirected
  invalid_gene_id <<- "2147483647"
  
  
  gmtOut <<- NA
  gmtCounter <<- 0
  gmtPathway <<- "dummy"
  
  ##count
  lowSizeMis <<- 0
  correctSizeMis <<- 0
  bigSizeMis <<- 0
  allLowSizeMis <<- data.frame(id=c(), pathway=c(), genes=c(), stringsAsFactors=F)
  allCorrectSizeMis <<- data.frame(id=c(), pathway=c(), genes=c(), stringsAsFactors=F)
  allBigSizeMis <<- data.frame(id=c(), pathway=c(), genes=c(), stringsAsFactors=F)
  allGis <<- data.frame(id=c(), pathway=c(), genes=c(), stringsAsFactors=F) 
  allGisGraphs <<- list() # giss are ordered (indexed) according to gisid in allGis
}

extractGisFromKgml <- function (kgmlNodeFile, kgmlEdgeFile){
  nodeData <- read.table(kgmlNodeFile, sep="\t", col.names=c("id","genes"), colClasses=c("numeric","character"))
  edgeData <- read.table(kgmlEdgeFile, col.names=c("from","to"), colClasses = c("numeric","numeric"))
  
  if(nrow(edgeData) <=0 )
    return (-1) 
  
  ## kegg network may have multiple edges between same pair of nodes.
  #print(kgmlNodeFile)
  print("removing multiple edges")
  print(Sys.time())
  edgeData <- removeMultipleEdges(edgeData, directed)
  print(Sys.time())
  print("removed")
  
  kgmlGraph <- buildIGraph(nodeData, edgeData, directed)
  components <- clusters(graph=kgmlGraph)
  
  # for each component, extract GIS
  for(comp in 1:components$no){
    
    componentNodes <- which(components$membership==comp)
    componentGraph <- induced.subgraph(graph=kgmlGraph, vids=componentNodes)
    
    componentNodeData <- nodeData[componentNodes, ,drop=F]
    componentNodeData$id <- 1:nrow(componentNodeData)
    componentEdgeData <- get.edgelist(graph=componentGraph, names=F)
    
    ##count low size network
    if(nrow(componentNodeData) < minNodeInGis){
      lowSizeMis <<- lowSizeMis + 1
    }
    else if(nrow(componentNodeData) <= maxNodeInGis){
      correctSizeMis <<- correctSizeMis + 1
    }
    
    
    gisLines <- extractGisFromConnectedGraph(componentNodeData, componentEdgeData)
    saveGisLines(gisLines)
  }
  
}
#debug(extractGisFromKgml)

removeMultipleEdges_v1 <- function(edgeData, directed=FALSE){
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

removeMultipleEdges_v2 <- function(edgeData, directed=FALSE){
  #print("debug-1")
  freq <- table(edgeData)
  rows <- rownames(freq)
  cols <- colnames(freq)
  print("debug-2")
  
  if(!directed){
    twoWayGenes <- intersect(rows, cols)
    if(length(twoWayGenes)>1){
      #print(length(twoWayGenes))
      print("debug-2.1")
      combs <- combn(twoWayGenes,2)
      print("debug-2.2")
      combs <- t(combs)
      print("debug-2.3")
      combinedFreq <- apply(combs, 1, function(pair){
          return(freq[pair[1],pair[2]] + freq[pair[2],pair[1]])
      })
      print("debug-2.4")
      freq[combs] <- combinedFreq
      print("debug-2.5")
      freq[combs[,c(2,1),drop=F]] <- 0
      
    }
    
  }
  
  print("debug-3")
  occured <- which(freq>0, arr.ind=T)
  relations <- apply(occured, 1, function(r){
    return(c(rows[r[1]], cols[r[2]]))
  })
  
  #print("debug-4")
  edgeData <- edgeData[1:nrow(occured),,drop=F]
  edgeData[,1] <- as.numeric(relations[1,])
  edgeData[,2] <- as.numeric(relations[2,])
  #print("debug-5")
  return(edgeData)
}

removeMultipleEdges <- function(edgeData, directed=FALSE){
  if(!directed){
    from <- apply(edgeData, 1, min)
    to <- apply(edgeData, 1, max)
    edgeData <- data.frame(from=from, to=to)
  }
  
  edgeData$from <- as.integer(edgeData$from)
  edgeData$to <- as.integer(edgeData$to)
  edgeData <- unique(edgeData)
  return(edgeData)
}

extractGisFromConnectedGraph <- function(nodeData, edgeData){
  if(nrow(nodeData)==0){
    #return("")
    miss = list()
    miss[[1]] = list(genes="", graph=graph.empty())
    return(miss)
  }
  
  kgmlGraph <- buildIGraph(nodeData, edgeData, directed)
  noOfNodes <- vcount(kgmlGraph)
  
  if(nrow(nodeData)<minNodeInGis){
    # taking gis from network having less than min nodes
    #return(genesInGraph(nodeData))
    miss = list()
    miss[[1]] = list(genes=genesInGraph(nodeData), graph=kgmlGraph)
    return(miss)
  }
  
  if(nrow(nodeData)<=maxNodeInGis){
    # the network is in range, take the full network
    #return(genesInGraph(nodeData))
    miss = list()
    miss[[1]] = list(genes=genesInGraph(nodeData), graph=kgmlGraph)
    return(miss)
  }
  
  
  ########## build the community dendrogram (graphList) ##########
  bet <- fastgreedy.community(kgmlGraph)
  
  graphList <- apply(nodeData, 1, function(node){
    g <- graph.empty(n=0, directed=F)
    g <- set.graph.attribute(graph=g, name="genes", value=node["genes"])
    g <- set.graph.attribute(graph=g, name="count", value=1)        # 1 node, count=1
    g <- set.graph.attribute(graph=g, name="nodes", value=as.integer(node["id"]))
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
  
  ######### find appropriate community ##########
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
  
  # last row in 'merges' is the first split, so loop from last to first row.
  lapply(nrow(bet$merges):1, findCommunityGraphs)
  
  
  # update discardedNodeList by adding discarded single nodes [graphIndex: 1:noOfNodes]
  takeDiscardedSingleNodes <- function(nodeIndex){
    if(!isTaken[nodeIndex]){
      discardedGraphs[[discardedGraphsCount+1]] <<- graphList[[nodeIndex]]
      discardedGraphsCount <<- discardedGraphsCount+1
    }
    return(NA)
  }
  lapply(1:noOfNodes, takeDiscardedSingleNodes)
  
  
  
  ############## merge discarded nodes with the closest community ############
  
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
  
  misGraphFromTakenGraphVariable <- function(tg){
    nodeIndexes <- getIntegerNodeVector(tg)
    nodeIndexes <- sort(x=nodeIndexes, decreasing=F)
    tgGraph <- induced.subgraph(graph=kgmlGraph, vids=nodeIndexes)
    tgNodeData <- nodeData[nodeIndexes, ,drop=F]
    tgNodeData$id <- 1:nrow(tgNodeData)
    tgGraph <- set.vertex.attribute(graph = tgGraph, name = "name", index = 1:nrow(tgNodeData), value = tgNodeData[, 2])
    return(tgGraph)
  }
  
  append.list <- function(l1, l2){
    nl <- l1
    for(i in 1:length(l2)){
      nl[[length(nl)+1]] <- l2[[i]]
    }
    return(nl)
  }
  
  # if nodes are discarded, assign them to some taken graph.
  if(length(discardedGraphs) > 0){
    takenNodes <- getAllNodes(takenGraphs)
    discardedNodes <- getAllNodes(discardedGraphs)
    
    while(length(discardedNodes) > 0){
      #print(paste(format(Sys.time(), "%X"), "#discarded nodes:",length(discardedNodes)))
      memberships <- sapply(discardedNodes, assignNode, takenNodes, 1)
      
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
    #print(paste(format(Sys.time(), "%X"), "#discarded nodes:",length(discardedNodes)))
    
  }
  
  
  ################ create MIS if size is within range, otherwise recursively generate MIS ##############
  
  if(takenGraphsCount==1){
    # could not divide more, so return it. Or divide in special way
    
    ##count
    bigSizeMis <<- bigSizeMis + 1
    
    #gisGenes <- genesInAnnotatedGraph(takenGraphs[[1]])
    #return(gisGenes)
    
    gline <- genesInAnnotatedGraph(takenGraphs[[1]])
    misGraph <- misGraphFromTakenGraphVariable(takenGraphs[[1]])
    miss = list()
    miss[[1]] = list(genes=gline, graph=misGraph)
    return(miss)
  }
  else{
    
    #gisGenes <- c()
    miss <- list()
    lapply(takenGraphs, function(tg){
      tg.count <- as.integer(get.graph.attribute(graph=tg, name="count"))
      if(tg.count<=maxNodeInGis){
        gline <- genesInAnnotatedGraph(tg)
        #gisGenes <<- c(gisGenes, gline)
        
        ## get the mis graph
        misGraph <- misGraphFromTakenGraphVariable(tg)
        newMiss <- list()
        newMiss[[1]] <- list(genes=gline, graph=misGraph)
        miss <<- append.list(miss, newMiss)
        
        ## count
        correctSizeMis <<- correctSizeMis + 1
      }
      else{
        
        #cat(paste("split network - #nodes:", tg.count, "\n"))
        
        nodeIndexes <- getIntegerNodeVector(tg)
        nodeIndexes <- sort(x=nodeIndexes, decreasing=F)
        
        tgGraph <- induced.subgraph(graph=kgmlGraph, vids=nodeIndexes)
        tgNodeData <- nodeData[nodeIndexes, ,drop=F]
        tgNodeData$id <- 1:nrow(tgNodeData)
        tgEdgeData <- get.edgelist(graph=tgGraph, names=F)
        
        #tgGisGenes <- extractGisFromConnectedGraph(tgNodeData, tgEdgeData)
        #gisGenes <<- c(gisGenes, tgGisGenes)
        tgMiss <- extractGisFromConnectedGraph(tgNodeData, tgEdgeData)
        miss <<- append.list(miss, tgMiss)
      }
    })
    
    #return(gisGenes)
    return(miss)
  }
}





buildIGraph <- function(nodeData, edgeData, directed=F){
  noOfNodes <- nrow(nodeData)
  edgeMatrix <- as.matrix(edgeData)
  edgeVector <- as.vector(t(edgeMatrix))
  gr <- graph(edges=edgeVector, n=noOfNodes, directed=directed)
  gr <- set.vertex.attribute(graph = gr, name = "name", index = 1:nrow(nodeData), value = nodeData[, 2])
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

getValidSortedUniqueGenes <- function(g.genes){
  g.genes <- strsplit(x=g.genes, split=",")[[1]]
  g.genes <- sort(unique(g.genes))
  g.genes <- setdiff(g.genes, invalid_gene_id)
  
  if(length(g.genes) == 0){
    return("")
  }
  
  line <- ""
  lapply(g.genes, function(gid){
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


saveGisLines <- function(miss){
  
  lapply(miss, function(mis){
    line = mis$genes
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
      allGisGraphs[[gmtCounter+1]] <<- mis$graph
      gmtCounter <<- gmtCounter+1
    }
    
    return(NA)
  })
  
  return(NA)
}

printInFiles <- function(){
  
  printInGmtFile()
  printInJsonFile()
  
}

printInGmtFile <- function(){
  
  apply(allGis, 1, function(gmtRow){
    line <- paste(str_trim(gmtRow["id"]), gmtRow["pathway"], gmtRow["genes"], sep="\t")
    write(file=gmtOut, x=line)  
    return(NA)
  })
  
}


printInJsonFile <- function(){
  
  misNetworks <- lapply(allGisGraphs, function(misGraph){
    return(misToJsonNetwork(misGraph))
  })
  
  
  if(length(misNetworks)>0){
    misList <- createEmptyMisList()
    for(i in 1:length(misNetworks)){
      # misNumber starts from 0 (as gisid starts from 0 in allGis)
      misList <- misList.addMis(misList, misNumber = i-1, mis = misNetworks[[i]])
    }
    
    # write in file
    cat(getJson(misList), file = jsonOut)
  }
  
}


misToJsonNetwork <- function(g){
  ## create empty network
  net <- createEmptyNetwork(dataAttrNames = c('label','expression', 'color'), 
                                     dataAttrTypes = c('string', 'string', 'string'), 
                                     dataAttrDefaultValues = c(NA,NA,'orange'),
                                     edgeAttrNames = c('type','directed'),
                                     edgeAttrTypes = c('string','boolean'),
                                     edgeAttrDefaultValues = list(NA, FALSE)
  );
  
  ## add each node
  misNodes <- get.vertex.attribute(g, name = "name")
  for(mnode in misNodes){
    net <- network.addNode(net, mnode , label=mnode)
  }
  
  ## add each edge
  misEdges <- get.edgelist(graph = g, names = T)
  if(nrow(misEdges) > 0){
    for(i in 1:nrow(misEdges)){
      net <- network.addEdge(net, as.character(i), source=misEdges[i,1], target=misEdges[i,2])
    }
  }
  
  return(net)
}



#debug(extractGisFromConnectedGraph)

generate_mis <- function(network, range=NA, gmtOutFile=paste0("results/", network, ".gmt"), jsonOutFile=paste0("results/", network, ".json")){
  
  mis_init()
  
  if(!(network %in% allowedNetworks)){
    stop("Please give one of the allowed networks (kegg/string).")
  }
  
  if(length(range)==1 && is.na(range)){
    range <- switch(network, kegg=c(5,15), string=c(5,25), pathwayapi=c(5,15))
  }
  
  if(class(range) !=  "numeric" || length(range)!=2 || any(range<=0) ){
    stop("Range must be a positive numeric array of length 2.")
  }
  
  if(network=="kegg"){
    cat("Please wait. It may take several minutes. Intel core-i5, 3.3GHz processor with 12GB ram takes about a minute.\n")
  }
  else if(network=="string"){
    cat("Please wait, it may take several minutes. Intel core-i5, 3.3GHz processor with 12GB ram takes about 15 minutes.\n")
  }
  else{
    cat("Please wait, it may take several minutes. Intel core-i5, 3.3GHz processor with 12GB ram takes about 40 minutes.\n")
  }
  
  #### setttings ##########
  network <<- network
  minNodeInGis <<- min(range)
  maxNodeInGis <<- max(range)
  
  #pathwaysFile <- switch(network, kegg="data/networks/kegg/kegg.txt", string="data/networks/string/string.txt")
  pathwaysFile <- paste0("data/networks/", network,"/", network,".txt")
  #gmtOutFile <- paste0("results/", network, ".gmt")
  #gmtOutFile <- outputFile
  
  gmtOut <<- file(gmtOutFile, open="w+")
  jsonOut <<- file(jsonOutFile, open="w+")
  gmtCounter <<- 0
  gmtPathway <<- "dummy"
  
  #### generate MIS from each pathway file ########
  pathways <- read.table(pathwaysFile)
  apply(pathways, 1, function(pathwayName){
    kgmlNodeFile <- paste("data/networks/", network, "/node-", pathwayName, ".txt", sep="")
    kgmlEdgeFile <- paste("data/networks/", network, "/edge-", pathwayName, ".txt", sep="")
    gmtPathway <<- pathwayName
    
    #cat(paste("processing file - ",gmtPathway,"\n"))
    extractGisFromKgml(kgmlNodeFile, kgmlEdgeFile)
    
  })
  
  printInFiles()
  cat("MISs are saved in 2 files:\n")  
  cat(paste("File-1: ", getwd(),"/", gmtOutFile, "\n", sep=""))  
  cat(paste("File-2: ", getwd(),"/", jsonOutFile, "\n", sep=""))  
  close(gmtOut)
  close(jsonOut)
}

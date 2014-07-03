############################################################################################
#### This script defines some methods to create an R object                            #####
#### which can be used as a cytoscape network.                                         #####
#### sample network:  http://cytoscapeweb.cytoscape.org/documentation/elements         #####
############################################################################################

library(jsonlite)
library(RCurl)

createDataSchemaObj <- function(dataAttrNames, dataAttrTypes, dataAttrDefaultValues=NULL, edgeAttrNames, edgeAttrTypes, edgeAttrDefaultValues=NULL){
  #if(class(dataAttrNames) != 'character' || class(dataAttrTypes) != 'character'|| class(edgeAttrNames) != 'character'|| class(edgeAttrTypes) != 'character')
  #  stop('Node and edge attribute names and types must be string.')
  
  if(length(dataAttrNames) != length(dataAttrTypes))
    stop('Node attribute names and values must be of same length.')
  
  if(length(edgeAttrNames) != length(edgeAttrTypes))
    stop('Edge attribute names and values must be of same length.')
  
  if(!is.null(dataAttrDefaultValues) && length(dataAttrNames) != length(dataAttrDefaultValues))
    stop('Node attribute default values will be either null or a vector/list of same length as attribute names.')
  
  if(!is.null(edgeAttrDefaultValues) && length(edgeAttrNames) != length(edgeAttrDefaultValues))
    stop('Edge attribute default values will be either null or a vector/list of same length as attribute names.')
  
  ## create node attributes
  schemaNodes <- list()
  dataLen <- length(dataAttrNames)
  if(dataLen > 0)
    schemaNodes <- lapply(1:dataLen, function(idx){
      nodeObj <- NULL
      if(!is.null(dataAttrDefaultValues) && !is.na(dataAttrDefaultValues[[idx]])){
        nodeObj <- list(name=dataAttrNames[idx], type=dataAttrTypes[idx], defValue=dataAttrDefaultValues[[idx]])
      } else{
        nodeObj <- list(name=dataAttrNames[idx], type=dataAttrTypes[idx])
      }
      return(nodeObj)
    })
  
  
  ## create edge attributes
  schemaEdges <- list()
  edgeLen <- length(edgeAttrNames)
  if(edgeLen>0)
    schemaEdges <- lapply(1:edgeLen, function(idx){
      edgeObj <- NULL
      if(!is.null(edgeAttrDefaultValues) && !is.na(edgeAttrDefaultValues[[idx]])){
        edgeObj <- list(name=edgeAttrNames[idx], type=edgeAttrTypes[idx], defValue=edgeAttrDefaultValues[[idx]])
      } else{
        edgeObj <- list(name=edgeAttrNames[idx], type=edgeAttrTypes[idx])
      }
      return(edgeObj)
    })
  
  ## create schema
  dataSchema <- list()
  dataSchema[['nodes']] <- schemaNodes
  dataSchema[['edges']] <- schemaEdges
  
  return(dataSchema)
}

createEmptyNetwork <- function(dataAttrNames, dataAttrTypes, dataAttrDefaultValues=NULL, edgeAttrNames, edgeAttrTypes, edgeAttrDefaultValues=NULL){
  dataSchema <- createDataSchemaObj(dataAttrNames, dataAttrTypes, dataAttrDefaultValues, edgeAttrNames, edgeAttrTypes, edgeAttrDefaultValues);
  dataValues <- list(nodes=list(), edges=list())
  network <- list(dataSchema=dataSchema, data=dataValues)
  return(network)  
}

network.addNode <- function(network, nodeId, ...){
  d <- list(id=nodeId)
  d <- merge.list(d, list(...))
  len <- length(network[['data']][['nodes']])
  network[['data']][['nodes']][[len+1]] <- d
  return(network)
}

network.addEdge <- function(network, edgeId, ...){
  d <- list(id=edgeId)
  d <- merge.list(d, list(...))
  len <- length(network[['data']][['edges']])
  network[['data']][['edges']][[len+1]] <- d
  return(network)
}

createEmptyMisList <- function(){
  return(list())
}

misList.addMis <- function(misList, misNumber, mis){
  misNode <- list(misnumber=misNumber, graph=mis)
  misList[[length(misList)+1]] <- misNode
  return(misList)
}

getJson <- function(obj){
  toJSON(obj, auto_unbox = T)
}

misToNetwork <- function(g){
  # create empty network
  
  
  # add each node
  
  # add each edge
  
}

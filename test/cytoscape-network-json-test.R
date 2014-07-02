source('cytoscape-network-json.R')

#### network creation example #######
emptyNetwork <- createEmptyNetwork(dataAttrNames = c('label','expression', 'color'), 
                              dataAttrTypes = c('string', 'string', 'string'), 
                              dataAttrDefaultValues = c(NA,NA,'orange'),
                              edgeAttrNames = c('type','directed'),
                              edgeAttrTypes = c('string','boolean'),
                              edgeAttrDefaultValues = list(NA, FALSE)
);

network1 <- network.addNode(emptyNetwork, '1', label="G1", expression="underexpressed", color="green")
network1 <- network.addNode(network1, '2', label="G2", expression="overexpressed", color="red")
network1 <- network.addNode(network1, '3', label="G3")

network1 <- network.addEdge(network1, '1', source='1', target='2')
network1 <- network.addEdge(network1, '2', source='1', target='3')
network1 <- network.addEdge(network1, '3', source='2', target='3')


network2 <- network.addNode(emptyNetwork, '1', label="H1")
network2 <- network.addNode(network2, '2', label="H2", expression="overexpressed", color='red')

network2 <- network.addEdge(network2, '1', source='1', target='2')

misList <- createEmptyMisList()
misList <- misList.addMis(misList, 1, network1)
misList <- misList.addMis(misList, 2, network2)

jsonObj <- getJson(misList)
cat(prettify(jsonObj))

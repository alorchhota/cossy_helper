## settings
pathwayApiDataDir <- 'data/pathwaycsv'
clusterOneInputFile <- "results\\clusterOneInputNetwork-PathwayAPI.txt"
gmtFile <- 'results\\pathwayapi_clustersone.gmt'

pathwayApiRelationsFile <- paste0(pathwayApiDataDir, '/pathway_rs.csv')
geneMappingFile <- paste0(pathwayApiDataDir, '/gene_mapping.csv')


## read input file
geneRelations <- read.csv(file=pathwayApiRelationsFile, header=T, sep=";", quote="\"",  stringsAsFactors=F, comment.char="")
geneMappings <- read.csv(file=geneMappingFile, header=T, sep=";", quote="\"",  stringsAsFactors=F, comment.char="", row.names=1)

## check if any gene name contains a white space
if(length(grep("\\s", geneMappings[,'gene_name']))>0)
  stop("gene name cannot have any white space")



## prepare input file for cluster one
clusterOneInputRelations <- geneRelations[,2:3,drop=F]
clusterOneInputRelations['weight'] <- 1

## convert gene_id to gene_names in gene-relations
getGeneName <- function(geneid){
  return(geneMappings[geneid,])
}
clusterOneInputRelations[,1] = getGeneName(clusterOneInputRelations[,1])
clusterOneInputRelations[,2] = getGeneName(clusterOneInputRelations[,2])

## remove any relation if it does not have any name
validRelations <- !is.na(clusterOneInputRelations[,1]) & !is.na(clusterOneInputRelations[,2])
clusterOneInputRelations <- clusterOneInputRelations[validRelations,]


## write the input file
write.table(clusterOneInputRelations, file=clusterOneInputFile, quote=F, sep=" ", row.names=F, col.names=F)

## invoke clusterONE
cmd <- paste0("java -jar cluster_one-1.0.jar ", clusterOneInputFile)
outputs <- system(cmd, intern=T)
clusters <- outputs[6:length(outputs)]

## output clusters in .gmt format (first entry in a line is clusterIndex, 2nd entry is name)
cat(paste(0:(length(clusters)-1), 'clusterone', clusters, sep='\t', collapse='\n'), file=gmtFile)


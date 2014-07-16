##########################################################################
####### This script downloads human kegg pathways                  #######
####### The list of human pathways is first generated              #######          
####### by searching hsa at http://www.genome.jp/kegg/pathway.html #######
##########################################################################

library(KEGGgraph)
library(KEGG.db)

pathwayListFile <- "data/raw/kegg-2014-07-16/kegg-human-pathways.txt"
outDir <- "results/kegg-2014-07-16/"


## save kegg db info and download date in readme.txt file.
dbinfo <- KEGG_dbInfo()
dlTime <- format(Sys.time(), "%Y-%m-%d %H:%M:%S KST")
dbinfo <- rbind(dbinfo, data.frame(name="Download Time", value=dlTime))
write.table(dbinfo, file = paste(outDir, "readme.txt"), quote = F, sep = ":\t", row.names = F, col.names = F)

## read pathway ids
pathwayTable <- read.table(file = pathwayListFile, header = T, sep = "\t", quote = "", stringsAsFactors = F)
pathwayIds <- pathwayTable$Entry
pathwayIds <- sapply(pathwayIds, function(pid) strsplit(pid,"hsa")[[1]][2])

## download files
downloadAndSaveKGML <- function(pId){
  dlfile <- paste0(outDir, "hsa", pId, ".kgml")
  retrieveKGML(pId, organism="hsa", destfile=dlfile, method="internal", quiet=TRUE)  
}

lapply(pathwayIds, downloadAndSaveKGML)

cat("Output saved in ", outDir)

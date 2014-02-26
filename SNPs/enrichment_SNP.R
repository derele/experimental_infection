
if(!exists("GT")){
    source("SNPs/genotyping.R")
}

library(topGO)


#######################################################
## This dapc stuff did not identify anything exciting
quite.sep.genes <- gsub("\\..*", "",
                       rownames(GT[dapc1$pca.loadings[,1]>0.006,]))

MF.seperate <- TOGO.all.onto("MF", names(transcript.2.GO), 
                             quite.sep.genes , transcript.2.GO)

GenTable(MF.seperate[[1]], MF.seperate[[2]])


BP.seperate <- TOGO.all.onto("BP", names(transcript.2.GO), 
                             quite.sep.genes , transcript.2.GO)

GenTable(BP.seperate[[1]], BP.seperate[[2]])

#######################################################
## dn/ds as such identified also nothing STRANGE FIXME!!! 
MF.pos.selected <- TOGO.all.onto("MF", names(transcript.2.GO), 
                             pos.selected , transcript.2.GO)

GenTable(MF.pos.selected[[1]], MF.pos.selected[[2]])


BP.pos.selected <- TOGO.all.onto("BP", names(transcript.2.GO), 
                                 pos.selected , transcript.2.GO)

GenTable(BP.pos.selected[[1]], BP.pos.selected[[2]])


#############################################################




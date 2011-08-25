library(ShortRead)
library(Rsamtools)
library(reshape)
library(DESeq)

## Rsamtools
what <- c("rname")# , "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(what = what)

files <- list.files("/home/ele/Data/RNAseq/mapping/", "*.bam$")

countsList <- list()
for (fi in files) {
  bam <- scanBam(paste("/home/ele/Data/RNAseq/mapping/", fi, sep=""), param=param)
  counts <- table(bam[[1]]$rname)
  counts.frame <- as.data.frame(counts)
  countsList[[fi]] <- counts.frame
}

counts.frame.long <- melt(countsList, id.vars="Var1")
counts.frame.long$variable <- NULL

counts.frame.wide <- reshape(counts.frame.long, v.names="value",
                             idvar="Var1", direction="wide", timevar="L1")
names(counts.frame.wide) <- gsub("^value\\.", "", names(counts.frame.wide))
names(counts.frame.wide) <- gsub("\\.bam$", "", names(counts.frame.wide))
row.names(counts.frame.wide) <- counts.frame.wide$Var1
counts.frame.wide$Var1 <- NULL

sex.conds <- factor(ifelse(grepl("M_\\d", names(counts.frame.wide)), "male", "female" ))
eel.conds <- factor(ifelse(grepl("^AA", names(counts.frame.wide)), "Aa", "Aj" ))
pop.conds <- factor(ifelse(grepl("\\w\\w_R.", names(counts.frame.wide)), "EU", "TW" ))

cds <- newCountDataSet(counts.frame.wide, sex.conds)
cds <- estimateSizeFactors(cds)
cds <- estimateVarianceFunctions(cds)
res <- nbinomTest(cds, "male", "female")

vsd <- getVarianceStabilizedData( cds )
dists <- dist( t( vsd ) )
idists <- as.matrix(dists)

## check that all tech-reps are clustering, m/f are clustering, etc...
heatmap (idists , symm=TRUE, margins = c (7,7))


### Sum up technical replicates!!!
counts.frame.long$L1 <- gsub("_\\d.bam", "", counts.frame.long$L1)
counts.frame.long$L1 <- as.factor(counts.frame.long$L1)
counts.frame <- as.data.frame(tapply(counts.frame.long$value,
                                     list(counts.frame.long$Var1, counts.frame.long$L1),
                                     sum))

sex.conds <- factor(ifelse(grepl("M$", names(counts.frame)), "male", "female" ))
eel.conds <- factor(ifelse(grepl("^AA", names(counts.frame)), "Aa", "Aj" ))
pop.conds <- factor(ifelse(grepl("\\w\\w_R.*", names(counts.frame)), "EU", "TW" ))


## Male/Female

cds <- newCountDataSet(counts.frame, sex.conds)
cds <- estimateSizeFactors(cds)
cds <- estimateVarianceFunctions(cds)
res <- nbinomTest(cds, "male", "female")

vsd <- getVarianceStabilizedData( cds )
dists <- dist( t( vsd ) )
idists <- as.matrix(dists)

## check that m/f are clustering, etc...
heatmap (idists , symm=TRUE, margins = c (7,7))

scvPlot(cds)

residualsEcdfPlot( cds, "male" )
residualsEcdfPlot( cds, "female" )

plot( 
     res$baseMean, 
     res$log2FoldChange, 
     log="x", pch=20, cex=.4, 
     col = ifelse( res$padj < .1, "red", "black" ),
     main="My data")

## sig <- res[res$padj < .01 & !is.na(res$padj),"id"]

## Eel

cds <- newCountDataSet(counts.frame, eel.conds)
cds <- estimateSizeFactors(cds)
cds <- estimateVarianceFunctions(cds)
res <- nbinomTest(cds, "Aa", "Aj")

scvPlot(cds)

residualsEcdfPlot( cds, "Aa" )
residualsEcdfPlot( cds, "Aj" )

plot( 
     res$baseMean, 
     res$log2FoldChange, 
     log="x", pch=20, cex=.4, 
     col = ifelse( res$padj < .1, "red", "black" ),
     main="My data")



## Populations

cds <- newCountDataSet(counts.frame, pop.conds)
cds <- estimateSizeFactors(cds)
cds <- estimateVarianceFunctions(cds)
res <- nbinomTest(cds, "EU", "TW")

scvPlot(cds)

residualsEcdfPlot( cds, "EU" )
residualsEcdfPlot( cds, "TW" )

plot( 
     res$baseMean, 
     res$log2FoldChange, 
     log="x", pch=20, cex=.4, 
     col = ifelse( res$padj < .1, "red", "black" ),
     main="My data")

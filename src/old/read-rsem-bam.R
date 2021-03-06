
## how many reads?
## command <- "ls /home/ele/Data/RNAseq/ | grep --perl \"AA|AJ_\" | parallel \'echo {}; grep \"@\" /home/ele/Data/RNAseq/{}/*_1.solfastq.gz.fastq | wc -l\'"
## nReads.raw <- readLines(pipe(command))
## nReads <- as.data.frame(t(matrix(nReads.raw, nrow=2)))
## names(nReads) <- c("library", "raw.reads")

## Rsamtools
library (Rsamtools)
what <- c("rname" , "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(what = what)

files <- list.files("/data/RNAseq/rsem_trinity/", "sorted.bam$")
## alternative method using library(GenomicRanges) would be:
## aligns <- readBamGappedAlignments("/home/ele/Data/RNAseq/mapping/AJ_T19M_1.bam")
## ie. bam <- scanBam("/home/ele/Data/RNAseq/mapping/AJ_T19M_1.bam", param=param)

countsList <- list()
for (fi in files) {
  bam <- scanBam(paste("/data/RNAseq/rsem_trinity/", fi, sep=""), param=param)
  counts <- table(bam[[1]]$rname)
  counts.frame <- as.data.frame(counts)
  countsList[[fi]] <- counts.frame
}

counts.frame.long <- melt(countsList, id.vars="Var1")
counts.frame.long$variable <- NULL

counts.frame.wide <- reshape(counts.frame.long, v.names="value",
                             idvar="Var1", direction="wide", timevar="L1")
names(counts.frame.wide) <- gsub("^value\\.", "", names(counts.frame.wide))
names(counts.frame.wide) <- gsub("\\.transcript\\.sorted\\.bam$", "", names(counts.frame.wide))
row.names(counts.frame.wide) <- counts.frame.wide$Var1
counts.frame.wide$Var1 <- NULL

## Just to see the clustering fo technical replicates!!!

sex.conds.raw <- factor(ifelse(grepl("\\dM$", names(counts.frame.wide)), "male", "female" ))
eel.conds.raw <- factor(ifelse(grepl("^AA", names(counts.frame.wide)), "AA", "AJ" ))
pop.conds.raw <- factor(ifelse(grepl("\\w\\w_R.*", names(counts.frame.wide)), "EU", "TW" ))
design.raw <- model.matrix(~sex.conds.raw*eel.conds.raw*pop.conds.raw)

cds.raw <- newCountDataSet(counts.frame.wide, design.raw)
cds.raw <- estimateSizeFactors(cds.raw)
cds.raw <- estimateDispersions(cds.raw)
vsd.raw <- getVarianceStabilizedData( cds.raw )
dists.raw <- dist( t( vsd.raw ) )
idists.raw <- as.matrix(dists.raw)

### Sum up technical replicates!!!
counts.frame.long$L1 <- gsub("_\\d.bam", "", counts.frame.long$L1)
counts.frame.long$L1 <- as.factor(counts.frame.long$L1)
counts.frame <- as.data.frame(tapply(counts.frame.long$value,
                                     list(counts.frame.long$Var1, counts.frame.long$L1),
                                     sum))

mapped.raw <- colSums(counts.frame)
## use counts frame for all calculations
### limit to good quality for the moment:
Ac.contigs <- contig.df[contig.df$Ac, "contig"]
counts.frame <- subset(counts.frame,
                       rownames(counts.frame)%in%Ac.contigs)
mapped.tax <- colSums(counts.frame)

MN.contigs <- contig.df[contig.df$AcMN, "contig"]
co.f <- subset(counts.frame,
               rownames(counts.frame)%in%MN.contigs &
               rowSums(counts.frame)>32)
mapped.good <- colSums(co.f)

map.tab <- merge(nReads, mapped.raw, by.x="library", by.y="row.names")
names(map.tab)[ncol(map.tab)] <- "raw.mapped"

map.tab <- merge(map.tab, mapped.tax, by.x="library", by.y="row.names")
names(map.tab)[ncol(map.tab)] <- "tax.mapped"

map.tab <- merge(map.tab, mapped.good, by.x="library", by.y="row.names")
names(map.tab)[ncol(map.tab)] <- "screened"

sex.conds <- factor(ifelse(grepl("M$", names(counts.frame)), "male", "female" ))
eel.conds <- factor(ifelse(grepl("^AA", names(counts.frame)), "Aa", "Aj" ))
pop.conds <- factor(ifelse(grepl("\\w\\w_R.*", names(counts.frame)), "EU", "TW" ))

design <- model.matrix(~sex.conds*eel.conds*pop.conds)



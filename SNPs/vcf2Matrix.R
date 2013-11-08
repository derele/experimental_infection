if(!exists("T.e")){
  source("/home/ele/thesis/experimental_infection/exp_test/edgeR.R")
}

VCF <- read.delim("/data/RNAseq/genotypes_above30.vcf",
                  comment.char="#",
                  header=FALSE,
                  as.is=TRUE)


## get rid of the non-Nematode tanscripts and of the low covered
## transcripts from which we do also not use expression evidence
VCF <- VCF[VCF$V1%in%rownames(T.e),]

rawGT <- VCF[,10:33]
rownames(rawGT) <- make.names(make.unique(VCF$V1))

## only one transcript per gene
rawGT <- rawGT[grepl("_seq1$", VCF$V1), ]


GT.per.gene.table <- table(gsub("\\.\\d+$", "", rownames(rawGT)))

names(rawGT) <- c("AA_R11M", "AA_R16M", "AA_R18F", "AA_R28F",
                  "AA_R2M", "AA_R8F", "AA_T12F", "AA_T20F", "AA_T24M", "AA_T3M",
                  "AA_T42M", "AA_T45F", "AJ_R1F", "AJ_R1M", "AJ_R3F", "AJ_R3M",
                  "AJ_R5F", "AJ_R5M", "AJ_T19M", "AJ_T20M", "AJ_T25M", "AJ_T26F",
                  "AJ_T5F", "AJ_T8F")

assign.GT <- function(x){
  ifelse(grepl("0/0", x), 0,
         ifelse(grepl("0/1", x), 1,
                ifelse(grepl("1/1", x), 2, NA)))
}

GT <- apply(rawGT, 2, assign.GT)

## hierarchical clustering
d <- dist(t(GT))
h <- hclust(d)
plot(h)
dev.off()


library(ape)
nj.GT <- nj(d)
plot(nj.GT, type="fan")
dev.off()

## parsimony tree
library(phangorn)
OC <- phyDat(as.matrix(t(GT)), type="USER", levels=c(0,1,2))
## OC <- acgt2ry(OC)
### dm.OC <- dist.logDet(OC)
tree.OC <- pratchet(OC)
plot(tree.OC)
dev.off()





## classical PCA
## pcrGT <- prcomp(t(GT))
## plot(pcrGT$loadings)
## dev.off()

## There could be an overabundance of alternative allele homozygotes
## in Europe
summary.factor(GT[,pop.conds%in%"TW"])
summary.factor(GT[,pop.conds%in%"EU"])

## Testing only for males/females shows that it is not detectable in
## the properly genotyped single-individual samples
summary.factor(GT[,pop.conds%in%"TW"&sex.conds%in%"female"])
summary.factor(GT[,pop.conds%in%"EU"&sex.conds%in%"female"])

summary.factor(GT[,pop.conds%in%"TW"&sex.conds%in%"male"])
summary.factor(GT[,pop.conds%in%"EU"&sex.conds%in%"male"])

per.sample.GT.sum <- t(apply(GT, 2, summary.factor))



## Some HardyWeinberg Stuff???
## apply(GT[2:3,], 1, function (x) table(x))

library(reshape)
library(adegenet)

## vignette("adegenet-genomics", package='adegenet')

GT.list <- apply(GT, 2, as.vector)
GT.ade <- new("genlight", GT.list)




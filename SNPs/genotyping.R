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

## only one transcript per gene
VCF <- VCF[grepl("_seq1$", VCF$V1), ]

rownames(VCF) <- make.names(make.unique(VCF$V1))

rawGT <- VCF[,10:33]

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
## give unique SNP-names
rownames(GT) <- rownames(rawGT)
## add .0 to the first SNP in a gene
rownames(GT) <- ifelse(grepl("\\.\\d+", rownames(GT)),
                       rownames(GT),
                       paste(rownames(GT), ".0", sep=""))

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
library(RSvgDevice)

## PCA would be nice
## vignette("adegenet-genomics", package='adegenet')
GT.list <- apply(GT, 1, as.vector)
GT.ade <- new("genlight", GT.list, ploidy=2)
GT.ade@pop <- pop.conds

### PCA
pca1 <- glPca(GT.ade, nf=4, n.cores=6)

devSVG("figures/geno_pca_col.svg")
myCol <- colorplot(pca1$scores, pca1$scores, transp=FALSE, cex=4)
abline(h=0,v=0, col="grey")
text(pca1$scores[,1], pca1$scores[,2], colnames(GT))
add.scatter.eig(pca1$eig, 2, 1, 2, posi="topright", inset=.05, ratio=.3)
dev.off()

## hierarchical clustering
d <- dist(t(GT))
h <- hclust(d)
## devSVG("figures/geno_hira.svg")
## plot(h)
## dev.off()

library(ape)
nj.GT <- nj(d)

devSVG("figures/geno_nj_col.svg")
plot(nj.GT, type="unrooted")
tiplabels(pch=20, col=myCol, cex=4)
dev.off()

## parsimony tree
library(phangorn)
OC <- phyDat(as.matrix(t(GT)), type="USER", levels=c(0,1,2))
## OC <- acgt2ry(OC)
### dm.OC <- dist.logDet(OC)
tree.OC <- pratchet(OC)

devSVG("figures/geno_pars_col.svg")
plot(tree.OC, type="unrooted")
tiplabels(pch=20, col=myCol, cex=4)
dev.off()

### Discriminant Analysis of Principal Components (DAPC)

clust.ade <- find.clusters(GT.ade)
dapc1 <- dapc(GT.ade, n.pca=5, n.da=1)

devSVG("figures/geno_dapc_discr.svg")
scatter(dapc1, scree.da=FALSE,
        bg="white",
        posi.pca="topright",
        legend=TRUE,
        txt.leg=levels(dapc1$grp),
        col=c("red","blue"))
dev.off()

## compoplot(dapc1,
##           col=c("red","blue"),
##           txt.leg=levels(dapc1$grp),
##           lab=colnames(GT),
##           posi="topright",
##           ncol=2)
## dev.off()

## loadingplot(dapc1$var.contr,
##             thres=1e-20)
## dev.off()

## seems that templock could be used for this...
temp <- seploc(GT.ade, block.size=1000, n.cores=6)


summary.factor(dapc1$var.contr<1e-15)
summary.factor(dapc1$pca.loadings[,1]==0)
## there are 203 SNPs with a very small variance contribution and a
## pca loading of 0 for the first axis. This identifies the most
## shared (less discriminating) SNPs across populations
## This logic somehow does not work!!!

pheatmap(GT[dapc1$pca.loadings[,1]>0.006,])
dev.off()

library(ggplot2)
GT[rowSums(GT[,pop.conds%in%"EU"])<1 & rowSums(GT[, pop.conds%in%"TW"])>15, ]
GT[rowSums(GT[,pop.conds%in%"EU"])>15 & rowSums(GT[, pop.conds%in%"TW"])<1, ]

amova(, distances, structures)

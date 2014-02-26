library(ggplot2)
library(reshape)
library(adegenet)
library(RSvgDevice)
library(phangorn)
library(ape)
library(RSvgDevice)

if(!exists("T.e")){
    source("exp_test/edgeR.R")
}

VCF <- read.delim("/data/A_crassus/RNAseq/genotypes_above30.vcf",
                  comment.char="#",
                  header=FALSE,
                  as.is=TRUE)

## get rid of the non-Nematode tanscripts and of the low covered
## transcripts from which we do also not use expression evidence
VCF <- VCF[VCF$V1%in%rownames(T.e),]

## only one transcript per gene
VCF <- VCF[grepl("_seq1$", VCF$V1), ]

rownames(VCF) <- paste(VCF$V1, ".", VCF$V2, ":",
                       VCF$V4, "-", VCF$V5, sep = "" )

rawGT <- VCF[,10:33]

GT.per.gene.table <- table(gsub("\\.\\d+.*", "", rownames(rawGT)))

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

## to avoid restriction of loci to genes later in making the as.loci
## or standard genind object
rownames(GT) <- gsub("\\.", "_" , rownames(GT))


n.worms.preped <- c(AA_R11M = 14, AA_R16M = 4, AA_R18F = 1, AA_R28F = 1,
                    AA_R2M = 4, AA_R8F = 1,  AA_T12F = 1, AA_T20F = 1,
                    AA_T24M = 3, AA_T3M = 4, AA_T42M = 1, AA_T45F = 1, 
                    AJ_R1F = 1, AJ_R1M = 1, AJ_R3F = 1, AJ_R3M = 2,
                    AJ_R5F = 1, AJ_R5M = 1, AJ_T19M = 7, AJ_T20M = 8,
                    AJ_T25M = 5, AJ_T26F = 1, AJ_T5F = 1, AJ_T8F = 1)

## Testing only for males/females shows that it is not detectable in
## the properly genotyped single-individual samples
table(GT[,pop.conds%in%"TW"&n.worms.preped==1])
table(GT[,pop.conds%in%"EU"&n.worms.preped==1])

table(GT[,pop.conds%in%"TW"&n.worms.preped>1])
table(GT[,pop.conds%in%"EU"&n.worms.preped>1])

per.sample.GT.sum <- as.data.frame(t(apply(GT, 2, table)))

per.sample.GT.sum$het.hom <-
    per.sample.GT.sum[, 2] /
    (per.sample.GT.sum[, 1] + per.sample.GT.sum[, 3]) 

per.sample.GT.sum <- cbind(per.sample.GT.sum, n.worms.preped)


## PCA would be nice
## vignette("adegenet-genomics", package='adegenet')
GT.list <- apply(GT, 1, as.vector)
GT.ade <- new("genlight", GT.list, ploidy=2,
              chromosome = gsub("\\.\\d+.*", "", rownames(GT)))
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

nj.GT <- nj(d)

devSVG("figures/geno_nj_col.svg")
plot(nj.GT, type="unrooted")
tiplabels(pch=20, col=myCol, cex=4)
dev.off()

## parsimony tree
OC <- phyDat(as.matrix(t(GT)), type="USER", levels=c(0,1,2))
## OC <- acgt2ry(OC)
### dm.OC <- dist.logDet(OC)
tree.OC <- pratchet(OC)

devSVG("figures/geno_pars_col.svg")
plot(tree.OC, type="unrooted")
tiplabels(pch=20, col=myCol, cex=4)
dev.off()

## cluster finding 
clust.ade <- find.clusters(GT.ade, n.clust = 2, n.pca = 5)

clust.4.ade <- find.clusters(GT.ade, n.clust = 4, n.pca = 7)

clust.opt.ade <- find.clusters(GT.ade, stat = "BIC", n.pca = 7)

dapc1 <- dapc(GT.ade, n.pca=5, n.da=1)

devSVG("figures/geno_dapc_discr.svg")
scatter(dapc1, scree.da=FALSE,
        bg="white",
        posi.pca="topright",
        legend=TRUE,
        txt.leg=levels(dapc1$grp),
        col=c("red","blue"))
dev.off()


## devSVG("figures/geno_dapc_col.svg")
## myCol <- colorplot(dapc1$scores, dapc1$scores, transp=FALSE, cex=4)
## abline(h=0,v=0, col="grey")
## text(pca1$scores[,1], pca1$scores[,2], colnames(GT))
## add.scatter.eig(pca1$eig, 2, 1, 2, posi="topright", inset=.05, ratio=.3)
## dev.off()

summary.factor(dapc1$var.contr<1e-15)
summary.factor(dapc1$pca.loadings[,1]==0)
## there are 203 SNPs with a very small variance contribution and a
## pca loading of 0 for the first axis. This identifies the most
## shared (less discriminating) SNPs across populations

## in the other direction it is possible to find SNPS that discrimitate perfectly
summary.factor(dapc1$pca.loadings[,1]>0.006)
summary.factor(dapc1$pca.loadings[,1]>0.01)

## This logic works!!! We can identify very few markers that seperate
## the populations well!! 
## pheatmap(GT[dapc1$pca.loadings[,1]>0.006,])
## dev.off()

## get per gene genotypes 
GT.gene <- as.data.frame(GT)
GT.gene$gene <- gsub("_seq\\d+\\.\\d+.*", "", rownames(GT))
GT.gene.list <- by(GT.gene, GT.gene$gene, function (x) {
    apply(x[, 1:24], 2, paste, collapse=".")})
GT.gene <- as.data.frame(do.call(rbind, GT.gene.list))

devSVG("figures/geno_most_224_heat.svg")
pheatmap(GT[dapc1$pca.loadings[,1]>0.01,])
dev.off()

GT.haplo <- GT>0

GT.haplo.diff <- sapply(0:11, function (i){ 
    rownames(GT.haplo[(rowSums(GT.haplo[, pop.conds%in%"EU"]) > i &
                       rowSums(GT.haplo[, pop.conds%in%"TW"]) == 0) |
                      (rowSums(GT.haplo[, pop.conds%in%"EU"]) == 0 &
                       rowSums(GT.haplo[, pop.conds%in%"TW"]) > i)
                      ,])})

DAPC.frame <- dapc1$pca.loadings
rownames(DAPC.frame) <- rownames(GT)

## Constructiong the full genind object for fst stats

GT.for.genind <- apply(GT, 2, function (x)
                       gsub("^2$", "2/2", gsub("^1$", "1/2", gsub("^0$", "1/1", x)))
                       )
rownames(GT.for.genind) <- rownames(GT)


## "expensive" Fst calculations 
if (!exists ("GT.ade.heavy") ) { 

    if( file.exists("/data/A_crassus/RNAseq/Fst.Rata")){
        load("/data/A_crassus/RNAseq/Fst.Rata")
    }

    else {
        GT.ade.heavy <- df2genind(t(GT.for.genind), ploidy=2, sep="/",
                                  pop=GT.ade@pop, ind.names=colnames(GT),
                                  loc.names=rownames(GT))

        ## Fst, Fis, Fit
        ## using hierfstat
        FST <- fstat(GT.ade.heavy[n.worms.preped==1, ])

        pair.matFst <- pairwise.fst(GT.ade.heavy[n.worms.preped==1, ],
                                    res.type="matrix", truenames=TRUE)

        ## use Fst from pegas
        GT.loci <- as.loci(GT.ade.heavy[n.worms.preped==1,])
        fsttab <- Fst(GT.loci)

        apply(fsttab, 2, mean, na.rm=TRUE)

        save(GT.ade.heavy, FST, pair.matFst, fsttab, GT.loci, 
             file = "/data/A_crassus/RNAseq/Fst.Rata")

        library(Rhh)
        het.table <- data.frame(ir(GT.ade.heavy@tab[n.worms.preped==1,]))
        rownames(het.table) <- GT.ade.heavy[n.worms.preped==1,]@ind.names

        het.table$hl <- hl(GT.ade.heavy@tab[n.worms.preped==1,])
        het.table$sh <- sh(GT.ade.heavy@tab[n.worms.preped==1,])

        ## expected heterozygosity again from adegenet
        GP.ade.heavy <- genind2genpop(GT.ade.heavy)
        exp.het <- Hs(GP.ade.heavy, truenames=TRUE)
    }
}



library(edgeR)
library(VennDiagram)
library(pheatmap)
library(RSvgDevice)
library(randomForest)

if(!exists("good.Tax.transripts")){
  source("src/2_taxon_screen.R")
}

## Reading the Rsem expression tables for genes
G.exp <- read.delim("/data/A_crassus/RNAseq/rsem_trinity/Trinity_genes.counts.matrix.TMM_normalized.FPKM")
rownames(G.exp) <- G.exp[,1]
G.exp$X <- NULL
names(G.exp) <- gsub("\\.genes\\.results", "", names(G.exp))
G.e <- round(G.exp)
G.e <- G.e[rowSums(G.e)>200, ]
G.e <- G.e[rownames(G.e)%in%good.Tax.genes,]
G.e <- as.matrix(G.e)

## Reading the Rsem expression tables for transcripts
T.exp <- read.delim("/data/A_crassus/RNAseq/rsem_trinity/Trinity_trans.counts.matrix.TMM_normalized.FPKM")
rownames(T.exp) <- T.exp[,1]
T.exp$X <- NULL
names(T.exp) <- gsub("\\.isoforms\\.results", "", names(T.exp))
T.e <- round(T.exp)
T.e <- T.e[rowSums(T.e)>100, ]
T.e <- T.e[rownames(T.e)%in%good.Tax.transripts,]

## to keep a record of what transcripts are in the analysis
# writeLines(rownames(T.e), "/data/A_crassus/RNAseq/Trinity/analysed_transcripts.txt")

## building the design matrix
sex.conds <- factor(ifelse(grepl("M$", names(G.exp)), "male", "female" ))
eel.conds <- factor(ifelse(grepl("^AA", names(G.exp)), "AA", "AJ" ))
pop.conds <- factor(ifelse(grepl("\\w\\w_R.*", names(G.exp)), "EU", "TW" ))

design <- model.matrix(~(sex.conds+eel.conds+pop.conds)^2)
rownames(design) <- colnames(T.e)

sample.annot <- as.data.frame(
    cbind(sex = as.character(sex.conds),
          host = as.character(eel.conds),
          worm.pop = as.character(pop.conds)))
rownames(sample.annot) <- colnames(T.e)

ed <- DGEList(G.e, lib.size=colSums(G.e))
td <- DGEList(T.e, lib.size=colSums(T.e))

devSVG("figures/mds_genes.svg")
plotMDS.DGEList(ed, top=500)
dev.off()

if(!file.exists("figures/heatmap_genes.pdf")){
    pdf("figures/heatmap_genes.pdf")
    pheatmap(ed$counts, scale="row", show_rownames = FALSE,
             annotation = sample.annot,
             )
    dev.off()
}

## exclude outlier samples
ed <- ed[, !rownames(ed$samples)%in%c("AJ_T26F", "AA_T42M")]
td <- td[, !rownames(td$samples)%in%c("AJ_T26F", "AA_T42M")]
design <- design[!rownames(design)%in%c("AJ_T26F", "AA_T42M"),]

## Data is normalized
## ed <- calcNormFactors(ed)
ed <- estimateGLMCommonDisp(ed, design=design)
ed <- estimateGLMTrendedDisp(ed, design=design)
ed <- estimateGLMTagwiseDisp(ed, design=design)

## Data is normalized
## td <- calcNormFactors(td)
td <- estimateGLMCommonDisp(td, design=design)
td <- estimateGLMTrendedDisp(td, design=design)
td <- estimateGLMTagwiseDisp(td, design=design)

g.glmfit<- glmFit(ed, design, dispersion=ed$tagwise.dispersion)
t.glmfit<- glmFit(td, design, dispersion=td$tagwise.dispersion)

glm.wrapper <- function(fit.obj, conds.regex){
  glm.list <- list()
  coe <- names(as.data.frame(fit.obj$design))
  coe.l <- lapply(conds.regex, function (x) grep(x, coe))
  for (i in 1:length(conds.regex)){
    glm.list[[conds.regex[[i]]]] <- glmLRT(fit.obj, coef=coe.l[[i]])
  }
  return(glm.list)
}

## selects all coefficients being contained in each other hierachically
combi.conds <- gsub(":", ".*", names(as.data.frame(g.glmfit$design))[2:7])

g.glm.l <- glm.wrapper(g.glmfit, combi.conds)
t.glm.l <- glm.wrapper(t.glmfit, combi.conds)

## topTags works (same as using p.adjust directly)
g.topTags.l <- lapply(g.glm.l, function (x){
  topTags(x, n=40000) ## set as high as the length 
})

t.topTags.l <- lapply(t.glm.l, function (x){
  topTags(x, n=40000) ## set as high as the length 
})

### Transcripts are the real thing; look much nicer!!!
t.test.l <- lapply(t.topTags.l, function (x)
                      rownames(x$table[x$table$FDR<0.05 &
                                       abs(x$table[,1]) > 1.5, ]))
lapply(t.test.l, length)

g.test.l <- lapply(g.topTags.l, function (x)
                      rownames(x$table[x$table$FDR<0.05 &
                                       abs(x$table[,1]) > 1.5, ]))
lapply(g.test.l, length)

devSVG("/home/ele/thesis/experimental_infection/figures/pop_heat_trans.svg")
pheatmap(td$counts[t.test.l[[3]], ], scale = "row" ,
         annotation = sample.annot, show_rownames = FALSE)
dev.off()

pdf("/home/ele/thesis/experimental_infection/figures/pop_heat_trans.pdf")
pheatmap(td$counts[t.test.l[[3]],], scale = "row" ,
         annotation = sample.annot, show_rownames = FALSE)
dev.off()

devSVG("/home/ele/thesis/experimental_infection/figures/eel_heat_trans.svg")
pheatmap(td$counts[t.test.l[[2]],], scale = "row" ,
         annotation = sample.annot, show_rownames = FALSE)
dev.off()

pdf("/home/ele/thesis/experimental_infection/figures/eel_heat_trans.pdf")
pheatmap(td$counts[t.test.l[[2]],], scale = "row" ,
         annotation = sample.annot, show_rownames = FALSE)
dev.off()


devSVG("figures/mds_trans_eel.svg")
plotMDS.DGEList(td[t.test.l[[2]],])
dev.off()

devSVG("figures/mds_trans_pop.svg")
plotMDS.DGEList(td[t.test.l[[3]],])
dev.off()


## using random forrest to see if contrasts can be destinguished 
library(randomForest)

rFo.sex <- randomForest(x=t(T.e),
                        sex.conds, subset = rownames(T.e)%in%colnames(td),
                        ntree=1000 )

rFo.eel <- randomForest(x=t(T.e),
                        y=eel.conds, subset = t.test.l[[2]],
                        ntree=1000)

rFo.pop <- randomForest(x=t(T.e),
                        y=pop.conds,
                        ntree=1000)

rFo.eel.pop <- randomForest(x=t(T.e),
                            y=eel.conds:pop.conds,
                            ntree=1000)

rFo.all <- randomForest(x=t(T.e), 
                        y=sex.conds:eel.conds:pop.conds,
                        ntree=1000)

sympa <- as.character(pop.conds:eel.conds)
sympa <- as.factor(ifelse(sympa%in%c("EU:AA", "TW:AJ"),
                          "sympa", "allopa"))

rFo.sympa <- randomForest(x=t(T.e),
                          y=sympa, subset = c(t.test.l[[3]], t.test.l[[2]]),
                          ntree=1000)

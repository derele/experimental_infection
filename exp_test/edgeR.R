library(edgeR)
library(VennDiagram)
library(pheatmap)

if(!exists ("b_annot")){
  source("/home/ele/thesis/experimental_infection/basic/trinity_annotation.R")
}

## Reading the Rsem expression tables for genes
G.exp <- read.delim("/data/RNAseq/rsem_trinity/Trinity_genes.counts.matrix.TMM_normalized.FPKM")
rownames(G.exp) <- G.exp[,1]
G.exp$X <- NULL
names(G.exp) <- gsub("\\.genes\\.results", "", names(G.exp))
G.e <- round(G.exp)
G.e <- G.e[rowSums(G.e)>200, ]

G.e <- G.e[rownames(G.e)%in%good.Tax.genes,]

## Reading the Rsem expression tables for transcripts
T.exp <- read.delim("/data/RNAseq/rsem_trinity/Trinity_trans.counts.matrix.TMM_normalized.FPKM")
rownames(T.exp) <- T.exp[,1]
T.exp$X <- NULL
names(T.exp) <- gsub("\\.isoforms\\.results", "", names(T.exp))
T.e <- round(T.exp)
T.e <- T.e[rowSums(T.e)>100, ]
T.e <- T.e[rownames(T.e)%in%good.Tax.transripts,]

## building the design matrix
sex.conds <- factor(ifelse(grepl("M$", names(G.exp)), "male", "female" ))
eel.conds <- factor(ifelse(grepl("^AA", names(G.exp)), "AA", "AJ" ))
pop.conds <- factor(ifelse(grepl("\\w\\w_R.*", names(G.exp)), "EU", "TW" ))

design <- model.matrix(~(sex.conds+eel.conds+pop.conds)^2)
design.wo.pop <- model.matrix(~(sex.conds+eel.conds)^2)
design.wo.eel <- model.matrix(~(sex.conds+pop.conds)^2)

ed <- DGEList(G.e, lib.size=colSums(G.e))
## Data is normalized
## ed <- calcNormFactors(ed)
ed <- estimateGLMCommonDisp(ed, design=design)
ed <- estimateGLMTrendedDisp(ed, design=design)
ed <- estimateGLMTagwiseDisp(ed, design=design)

td <- DGEList(T.e, lib.size=colSums(T.e))
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
t.contigs.l <- lapply(t.topTags.l, function (x)
                      rownames(x$table[x$table$FDR<0.2 &
                                       abs(x$table[,1]) > 0, ]))
lapply(t.contigs.l, length)

g.contigs.l <- lapply(g.topTags.l, function (x)
                      rownames(x$table[x$table$FDR<0.2 &
                                       abs(x$table[,1]) > 0, ]))
lapply(g.contigs.l, length)

pop.straight.t <- t.contigs.l[[3]][!t.contigs.l[[3]]%in%t.contigs.l[[6]]]

## heatmap.4.contigs <- function (pdata, ...){
##   ## pdata <- merge(pdata,
##   ##                b_annot[b_annot$seqnames%in%rownames(pdata),
##   ##                        c("seqnames", "nuc_nr_annot")],
##   ##                by.x = 0, by.y = "seqnames",
##   ##                all.x = TRUE)
##   ## rownames(pdata) <- make.unique(make.names(pdata$nuc_nr_annot))
##   ## pdata$Row.names <- NULL
##   ## pdata$nuc_nr_annot <- NULL
##   annot <- data.frame(pop.conds,
##                  sex.conds,
##                  eel.conds)
##   rownames(annot) <- colnames(T.e)
##   rownames(pdata) <- NULL
##   pheatmap(pdata, annotation = annot, scale="row", ...)
## }

## pdf("/home/ele/thesis/experimental_infection/figures/pop_straight_heat.pdf")
## heatmap.4.contigs(T.e[pop.straight.t,] )
## dev.off()

## pdf("/home/ele/thesis/experimental_infection/figures/pop_heat.pdf")
## heatmap.4.contigs(T.e[t.contigs.l[[3]],])
## dev.off()

## ## heatmap.4.contigs(G.e[g.contigs.l[[2]], ])

## ## use Edge DEseq to get variance stabilized data ??
## ### vsd <- getVarianceStabilizedData( cds.DE )

## ## eel.conds <- eel.conds[sex.conds%in%"male"]
## ## pop.conds <- pop.conds[sex.conds%in%"male"]
## PFC.full <- as.data.frame(predFC(T.e, # [,sex.conds%in%"male"],
##                                  design = design,
##                                  dispersion=td$tagwise.dispersion))
## colnames(PFC.full) <- gsub(":", ".", colnames(PFC.full))
## colnames(PFC.full) <- gsub(":", ".", colnames(PFC.full))
## colnames(PFC.full) <- gsub("\\(|\\)", "", colnames(PFC.full))

## ## reset gene-expression to normal values in the other eel-species
## EelVSEelPop <- ggplot(PFC.full, aes(eel.condsAJ + pop.condsTW,
##                             eel.condsAJ  +  eel.condsAJ.pop.condsTW,
##                                     color = rownames(PFC.full)%in%t.contigs.l[[6]])) +
##   geom_point(alpha = 0.1) +
##   geom_smooth(method = "lm", se=FALSE, color="darkgrey", alpha = 0.5) +
##   scale_x_continuous("log10 fold change European vs. Taiwanese host") +
##   scale_y_continuous("log10 fold European vs. Taiwanese host with interaction effect ") + 
##   scale_color_manual("significance\nof interaction", values = c("black", "red"))+
##   theme_bw() +
##   geom_abline(yintercept=0,slope=1)

## ## the differentially expressed genes for the pop*eel interaction
## ## reset gene-expression to normal values in the other eel-species
## ## all the population-different genes are in setting back
## ## eel-differences
## ## EelVSEelPop <- EelVSEelPop +
## ##   geom_point(data = PFC.full[t.contigs.l[[3]],],
## ##              aes(eel.condsAJ  + pop.condsTW ,
## ##                  eel.condsAJ + pop.condsTW + eel.condsAJ.pop.condsTW),
## ##              color = "indianred", size = 2) +
## ##   geom_smooth(data = PFC.full[t.contigs.l[[3]],], method = "lm", se=FALSE, color="indianred")

## ## However, the significant eel-differences are  reverted 
## library(RSvgDevice)
## devSVG(file="/home/ele/Dropbox/svg_work/main_vs_inter_working.svg")
## EelVSEelPop + 
##   geom_point(data = PFC.full[t.contigs.l[[6]],],
##              aes(eel.condsAJ + pop.condsTW,
##                  eel.condsAJ + eel.condsAJ.pop.condsTW), color="red", size = 2) + 
##   geom_smooth(data = PFC.full[t.contigs.l[[6]],], method = "lm", se=FALSE, color="red")
## dev.off()


## ## pheatmap(G.e[g.contigs.l[[3]], ], scale="row")

## ## contigs.all <- rownames(T.e)
## ## library(VennDiagram)

## ## foo <- list(sex_all = match(t.contigs.l[[1]], contigs.all),
## ##             eel_all = match(t.contigs.l[[2]], contigs.all),
## ##             pop_all = match(t.contigs.l[[3]], contigs.all)
## ##             )

## contigs.all <- rownames(T.e)

## bar <- as.matrix(data.frame(
##   sex_all = as.numeric(contigs.all%in%t.contigs.l[[1]]),
##   eel_all = as.numeric(contigs.all%in%t.contigs.l[[2]]),
##   pop_all = as.numeric(contigs.all%in%t.contigs.l[[3]])
##   ))

## library(limma)

## devSVG(file="/home/ele/Dropbox/svg_work/venn.svg")
## vennDiagram(bar)
## dev.off()

## ## venn.diagram(foo,
## ## filename = "test.tiff")

## ## ## male only

## ## Me <- G.e[, sex.conds%in%"male"]

## ## heatmap.4.contigs(Me[g.contigs.l[[2]], ])
## ## heatmap.4.contigs(Me[g.contigs.l[[3]], ])

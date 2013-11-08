
if(!exists ("T.e")){
  source("/home/ele/thesis/experimental_infection/basic/edgeR.R")
}


T.e.M <- T.e[,sex.conds%in%"male"]
T.e.F <- T.e[,sex.conds%in%"female"]

## building the design matrix
eel.conds.M <- factor(ifelse(grepl("^AA", names(T.e.M)), "AA", "AJ" ))
pop.conds.M <- factor(ifelse(grepl("\\w\\w_R.*", names(T.e.M)), "EU", "TW" ))

## building the design matrix
eel.conds.F <- factor(ifelse(grepl("^AA", names(T.e.F)), "AA", "AJ" ))
pop.conds.F <- factor(ifelse(grepl("\\w\\w_R.*", names(T.e.F)), "EU", "TW" ))

design.M <- model.matrix(~(eel.conds.M+pop.conds.M)^2)
design.F <- model.matrix(~(eel.conds.F+pop.conds.F)^2)

td.M <- DGEList(T.e.M, lib.size=colSums(T.e.M))
## Data is normalized
## td.M <- calcNormFactors(td.M)
td.M <- estimateGLMCommonDisp(td.M, design=design.M)
td.M <- estimateGLMTrendedDisp(td.M, design=design.M)
td.M <- estimateGLMTagwiseDisp(td.M, design=design.M)

td.F <- DGEList(T.e.F, lib.size=colSums(T.e.F))
## Data is normalized
## td.F <- calcNormFactors(td.F)
td.F <- estimateGLMCommonDisp(td.F, design=design.F)
td.F <- estimateGLMTrendedDisp(td.F, design=design.F)
td.F <- estimateGLMTagwiseDisp(td.F, design=design.F)

t.glmfit.M<- glmFit(td.M, design.M, dispersion=td.M$tagwise.dispersion)
t.glmfit.F<- glmFit(td.F, design.F, dispersion=td.F$tagwise.dispersion)

## selects all coefficients being contained in each other hierachically
combi.conds.M <- gsub(":", ".*", names(as.data.frame(t.glmfit.M$design))[2:4])
combi.conds.F <- gsub(":", ".*", names(as.data.frame(t.glmfit.F$design))[2:4])

t.glm.l.M <- glm.wrapper(t.glmfit.M, combi.conds.M)
t.glm.l.F <- glm.wrapper(t.glmfit.F, combi.conds.F)

t.topTags.l.M <- lapply(t.glm.l.M, function (x){
  topTags(x, n=40000) ## set as high as the length 
})

t.topTags.l.F <- lapply(t.glm.l.F, function (x){
  topTags(x, n=40000) ## set as high as the length 
})

### Transcripts are the real thing; look much nicer!!!
t.contigs.l.M <- lapply(t.topTags.l.M, function (x)
                        rownames(x$table[x$table$FDR<0.2 &
                                         abs(x$table[,1]) > 0, ]))

### Transcripts are the real thing; look much nicer!!!
t.contigs.l.F <- lapply(t.topTags.l.F, function (x)
                        rownames(x$table[x$table$FDR<0.2 &
                                         abs(x$table[,1]) > 0, ]))

bar <- as.matrix(data.frame(
  pop_all = as.numeric(contigs.all%in%t.contigs.l.M[[1]]),
  eel_all = as.numeric(contigs.all%in%t.contigs.l.M[[2]])
  ))

library(limma)
devSVG(file="/home/ele/Dropbox/svg_work/venn_M.svg")
vennDiagram(bar)
dev.off()

bar <- as.matrix(data.frame(
  pop_all = as.numeric(contigs.all%in%t.contigs.l.F[[1]]),
  eel_all = as.numeric(contigs.all%in%t.contigs.l.F[[2]])
  ))

library(limma)
devSVG(file="/home/ele/Dropbox/svg_work/venn_F.svg")
vennDiagram(bar)
dev.off()


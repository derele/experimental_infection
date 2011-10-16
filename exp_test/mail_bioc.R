## generate a df of neg. bionm. counts
y <- as.data.frame(matrix(rnbinom(6000,mu=10,size=10),ncol=24))
names(y) <- c("AA_R11M", "AA_R16M", "AA_R18F", "AA_R28F",
              "AA_R2M",  "AA_R8F",  "AA_T12F", "AA_T20F",
              "AA_T24M", "AA_T3M", "AA_T42M", "AA_T45F",
              "AJ_R1F",  "AJ_R1M",  "AJ_R3F",  "AJ_R3M",
              "AJ_R5F",  "AJ_R5M",  "AJ_T19M", "AJ_T20M",
              "AJ_T25M", "AJ_T26F", "AJ_T5F",  "AJ_T8F")

sex.conds <- factor(ifelse(grepl("M$", names(y)), "male", "female" ))
eel.conds <- factor(ifelse(grepl("^AA", names(y)), "Aa", "Aj" ))
pop.conds <- factor(ifelse(grepl("\\w\\w_R.*", names(y)), "EU", "TW" ))

design <- model.matrix(~sex.conds*eel.conds*pop.conds)

## Counts frame to full DGEList
d <- DGEList(y, lib.size=colSums(y))
d <- calcNormFactors(d)
d <- estimateGLMCommonDisp(d, design=design)
d <- estimateGLMTrendedDisp(d, design=design)
d <- estimateGLMTagwiseDisp(d, design=design)

glmfit <- glmFit(d, design, dispersion=d$tagwise.dispersion)

glm.wrapper <- function(de.obj, fit.obj, conds.regex){
  glm.list <- list()
  coe <- names(as.data.frame(fit.obj$design))
  coe.l <- lapply(conds.regex, function (x) grep(x, coe))
  for (i in 1:length(conds.regex)){
    glm.list[[conds.regex[[i]]]] <- glmLRT(de.obj, fit.obj, coef=coe.l[[i]])
  }
  return(glm.list)
}

## selects all coefficients being contained in each other hierachically
combi.conds <- gsub(":", ".*", names(as.data.frame(glmfit$design))[2:8])
glm.l <- glm.wrapper(d, glmfit, combi.conds)

## show what is compared
lapply(glm.l, function (x) x$comparison)

## topTags works (same as using p.adjust directly)
topTags.l <- lapply(glm.l, function (x){
  tt <- topTags(x, n=40000) ## set as high as the length 
  tt[tt$table$adj.P.Val<0.05] ## get only below adj.P
})

## decideTestsDGE does not work
decideTestsDGE.l <- lapply(glm.l, function (x){
   subset(x$table, (decideTestsDGE(x, p.value = .05))!=0)})


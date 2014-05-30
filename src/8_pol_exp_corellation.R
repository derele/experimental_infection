if(!exists("VAR")){
    source("src/5_coding_poly.R")
}

if(!exists("transcript.2.GO")){
    source("src/6_trinity_annotation.R")
}

if(!exists("contig.geno")){
    source("src/7_enrichment.R")
}


## 
contig.geno$Row.names <- NULL

## Correlation of differential expression and discrimination between 
contig.geno.sex <- merge(contig.geno, t.topTags.l[[1]]$table, by = 0)

## expression 
cor.test(contig.geno.sex$max.dapc, abs(contig.geno.sex$logFC.sex.condsmale),
         method="kendall")
cor.test(contig.geno.sex$max.dapc, contig.geno.sex$FDR, method="kendall")
cor.test(contig.geno.sex$max.dapc, contig.geno.sex$PValue, method="kendall")
## more different for higher geneotypic differenciation

cor.test(contig.geno.sex$Fst.max, abs(contig.geno.sex$logFC.sex.condsmale),
         method="kendall")
cor.test(contig.geno.sex$Fst.max, log(contig.geno.sex$FDR), method="kendall")
cor.test(contig.geno.sex$Fst.max, log(contig.geno.sex$PValue), method="kendall")

cor.test(contig.geno.sex$Fst.mean, abs(contig.geno.sex$logFC.sex.condsmale),
         method="kendall")
cor.test(contig.geno.sex$Fst.mean, log(contig.geno.sex$FDR), method="kendall")
cor.test(contig.geno.sex$Fst.mean, log(contig.geno.sex$PValue), method="kendall")

cor.test(contig.geno.sex$contig.dn.ds, abs(contig.geno.sex$logFC.sex.condsmale),
         method="kendall")
cor.test(contig.geno.sex$contig.dn.ds, contig.geno.sex$FDR,
         method="kendall")
cor.test(contig.geno.sex$contig.dn.ds, contig.geno.sex$PValue,
         method="kendall")

################### Eel  #############################################
contig.geno.eel <- merge(contig.geno, t.topTags.l[[2]]$table, by = 0)

cor.test(contig.geno.eel$max.dapc, abs(contig.geno.eel$logFC.eel.condsAJ),
         method="kendall")
cor.test(contig.geno.eel$max.dapc, log(contig.geno.eel$FDR), method="kendall")
cor.test(contig.geno.eel$max.dapc, log(contig.geno.eel$PValue), method="kendall")

cor.test(contig.geno.eel$Fst.max, abs(contig.geno.eel$logFC.eel.condsAJ),
         method="kendall")
cor.test(contig.geno.eel$Fst.max, log(contig.geno.eel$FDR), method="kendall")
cor.test(contig.geno.eel$Fst.max, log(contig.geno.eel$PValue), method="kendall")

cor.test(contig.geno.eel$Fst.mean, abs(contig.geno.eel$logFC.eel.condsAJ),
         method="kendall")
cor.test(contig.geno.eel$Fst.mean, log(contig.geno.eel$FDR), method="kendall")
cor.test(contig.geno.eel$Fst.mean, log(contig.geno.eel$PValue), method="kendall")

cor.test(contig.geno.eel$contig.dn.ds, abs(contig.geno.eel$logFC.eel.condsAJ),
         method="kendall")
cor.test(contig.geno.eel$contig.dn.ds, contig.geno.eel$FDR,
         method="kendall")
cor.test(contig.geno.eel$contig.dn.ds, contig.geno.eel$PValue,
         method="kendall")

################### Population #############################################
## dapc with all
contig.geno.pop <- merge(contig.geno, t.topTags.l[[3]]$table, by = 0)

cor.test(contig.geno.pop$max.dapc, abs(contig.geno.pop$logFC.pop.condsTW),
         method="kendall")
cor.test(contig.geno.pop$max.dapc, log(contig.geno.pop$FDR), method="kendall")
cor.test(contig.geno.pop$max.dapc, log(contig.geno.pop$PValue), method="kendall")

## Fst with all 
cor.test(contig.geno.pop$Fst.max, abs(contig.geno.pop$logFC.pop.condsTW),
         method="kendall")
cor.test(contig.geno.pop$Fst.max, log(contig.geno.pop$FDR), method="kendall")
cor.test(contig.geno.pop$Fst.max, log(contig.geno.pop$PValue), method="kendall")

cor.test(contig.geno.pop$Fst.mean, abs(contig.geno.pop$logFC.pop.condsTW),
         method="kendall")
cor.test(contig.geno.pop$Fst.mean, log(contig.geno.pop$FDR), method="kendall")
cor.test(contig.geno.pop$Fst.mean, log(contig.geno.pop$PValue), method="kendall")

## dnds with all 
cor.test(contig.geno.pop$contig.dn.ds, abs(contig.geno.pop$logFC.pop.condsTW),
         method="kendall")
cor.test(contig.geno.pop$contig.dn.ds, contig.geno.pop$FDR,
         method="kendall")
cor.test(contig.geno.pop$contig.dn.ds, contig.geno.pop$PValue,
         method="kendall")


## But what about the influece of overall expression

## logCPM is the overall expression: Ha! Not correlated
all(contig.geno.sex$logCPM == contig.geno.pop$logCPM)
all(contig.geno.sex$logCPM == contig.geno.eel$logCPM)

cor.test(contig.geno.sex$max.dapc, contig.geno.sex$logCPM,
         method="kendall")

cor.test(contig.geno.sex$mean.dapc, contig.geno.sex$logCPM,
         method="kendall")

cor.test(contig.geno.sex$Fst.max, contig.geno.sex$logCPM,
         method="kendall")

cor.test(contig.geno.sex$Fst.mean, contig.geno.sex$logCPM,
         method="kendall")

### and as a test FDR is correlated with logCPM as expected
cor.test(contig.geno.sex$FDR, contig.geno.sex$logCPM,
         method="kendall")

## The interesting correlation structure of dnds with overall
## exprssion and number of SNPs
cor.test(contig.geno.pop$contig.dn.ds, contig.geno.pop$logCPM,
         method="kendall")

cor.test(contig.geno.pop$n.SNPs, contig.geno.pop$logCPM,
         method="kendall")

cor.test(contig.geno.pop$n.SNPs, contig.geno.pop$contig.dn.ds,
         method="kendall")


dn.ds.mod <- glm(log(contig.dn.ds) ~n.SNPs+logCPM+log(FDR)+
                 Fst.mean+Fst.max+
                 logFC.sex.condsmale.pop.condsTW,
                 data=contig.geno.sex, subset=contig.dn.ds!=4 & contig.dn.ds>0)

## nemove n.SNPs
dn.ds.mod <- glm(log(contig.dn.ds) ~ logCPM + log(FDR)+
                 Fst.mean+Fst.max+
                 logFC.sex.condsmale.pop.condsTW,
                 data=contig.geno.sex, subset=contig.dn.ds!=4 & contig.dn.ds>0)

## remove logFC.sex.condsmale.pop.condsTW
dn.ds.mod <- glm(log(contig.dn.ds) ~ logCPM + log(FDR),
                 data=contig.geno.sex, subset=contig.dn.ds!=4 & contig.dn.ds>0)

summary(dn.ds.mod)

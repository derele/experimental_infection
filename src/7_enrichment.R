if(!exists("VAR")){
    source("src/5_coding_poly.R")
}

if(!exists("transcript.2.GO")){
    source("src/6_trinity_annotation.R")
}

library(topGO)
library(plyr)


MF.sex <- TOGO.all.onto("MF", names(transcript.2.GO), 
                        t.test.l[[1]], transcript.2.GO)
GenTable(MF.sex[[1]], MF.sex[[2]])

BP.sex <- TOGO.all.onto("BP", names(transcript.2.GO), 
                        t.test.l[[1]], transcript.2.GO)
GenTable(BP.sex[[1]], BP.sex[[2]])

CC.sex <- TOGO.all.onto("CC", names(transcript.2.GO), 
                         t.test.l[[1]], transcript.2.GO)
GenTable(CC.sex[[1]], CC.sex[[2]])
##############

MF.eel <- TOGO.all.onto("MF", names(transcript.2.GO), 
                        t.test.l[[2]], transcript.2.GO)
GenTable(MF.eel[[1]], MF.eel[[2]])

BP.eel <- TOGO.all.onto("BP", names(transcript.2.GO), 
                        t.test.l[[2]], transcript.2.GO)
GenTable(BP.eel[[1]], BP.eel[[2]])

CC.eel <- TOGO.all.onto("CC", names(transcript.2.GO), 
                         t.test.l[[2]], transcript.2.GO)
GenTable(CC.eel[[1]], CC.eel[[2]])
##############

MF.pop <- TOGO.all.onto("MF", names(transcript.2.GO), 
                        t.test.l[[3]], transcript.2.GO)
GenTable(MF.pop[[1]], MF.pop[[2]])

BP.pop <- TOGO.all.onto("BP", names(transcript.2.GO), 
                        t.test.l[[3]], transcript.2.GO)
GenTable(BP.pop[[1]], BP.pop[[2]])

CC.pop <- TOGO.all.onto("CC", names(transcript.2.GO), 
                         t.test.l[[3]], transcript.2.GO)
GenTable(CC.pop[[1]], CC.pop[[2]])


#######################################################
## This dapc stuff did not identify anything exciting
quite.sep.genes <-  rownames(GT[dapc1$var.contr[,1]>0.0001,])
quite.sep.genes <- unique(gsub("_\\d+:.*", "", quite.sep.genes))

MF.seperate <- TOGO.all.onto("MF", names(transcript.2.GO), 
                             quite.sep.genes , transcript.2.GO)
GenTable(MF.seperate[[1]], MF.seperate[[2]])

BP.seperate <- TOGO.all.onto("BP", names(transcript.2.GO), 
                        quite.sep.genes , transcript.2.GO)

GenTable(BP.seperate[[1]], BP.seperate[[2]])

CC.seperate <- TOGO.all.onto("CC", names(transcript.2.GO), 
                        quite.sep.genes , transcript.2.GO)

GenTable(CC.seperate[[1]], CC.seperate[[2]])

#######################################################
## dn/ds as such identified also nothing STRANGE FIXME!!! 
MF.pos.selected <- TOGO.all.onto("MF", names(transcript.2.GO), 
                                 pos.selected , transcript.2.GO)

GenTable(MF.pos.selected[[1]], MF.pos.selected[[2]])


BP.pos.selected <- TOGO.all.onto("BP", names(transcript.2.GO), 
                                 pos.selected , transcript.2.GO)

GenTable(BP.pos.selected[[1]], BP.pos.selected[[2]])


CC.pos.selected <- TOGO.all.onto("CC", names(transcript.2.GO), 
                                 pos.selected , transcript.2.GO)

GenTable(CC.pos.selected[[1]], CC.pos.selected[[2]])
#############################################################


pdf("figures/fst.dapc.pdf")
ggplot(VAR, aes(x=Fst, y=dapc.var)) + geom_point(alpha = 0.1) + geom_density2d()
dev.off()

pdf("figures/fst.ax1.pdf")
ggplot(VAR, aes(x=Fst, y=abs(Axis1))) + geom_point(alpha = 0.1) + geom_density2d()
dev.off()

pdf("figures/dapc.ax1.pdf")
ggplot(VAR, aes(x=dapc.var, y=abs(Axis1))) + geom_point(alpha = 0.1) + geom_density2d()
dev.off()


VARsum <- melt(VAR[, c("effect", "Fst", "Fis", "dapc.var", "Axis1")]) 

cast(VARsum, effect ~ variable , mean, na.rm = TRUE,
     subset=variable%in%c("Fst", "Fis", "dapc.var", "Axis1"))

devSVG("figures/Fst_Fis_dens.svg", width=14, height=7)
ggplot(subset(VARsum, (!is.na(effect) & !effect%in%"Nonsense") &
              variable%in%c("Fst", "Fis")),
       aes(value, color=effect)) +
    facet_wrap(~variable) +
    geom_density() +
    theme_bw() +
    scale_x_log10(breaks=1/10^(1:5))
dev.off()

devSVG("figures/dapc_axis_dens.svg", width=14, height=7)
ggplot(subset(VARsum, (!is.na(effect) & !effect%in%"Nonsense") &
              variable%in%c("dapc.var", "Axis1")),
       aes(value, color=effect)) +
    facet_wrap(~variable, scales="free") +
    geom_density() +
    theme_bw() +
    scale_x_log10()
dev.off()

contig.geno <- cbind(mean.dapc=tapply(VAR$dapc.var,
                         as.character(VAR$seqnames),
                         function (x) mean(abs(x))))

contig.geno <- cbind(contig.geno, max.dapc=tapply(VAR$dapc.var,
                                      as.character(VAR$seqnames),
                                      function (x) max(abs(x))))

contig.geno <- cbind(contig.geno,
                     Fst.mean=tapply(VAR$Fst,
                         as.character(VAR$seqnames), mean, na.rm=TRUE))

contig.geno <- cbind(contig.geno,
                     Fst.max=tapply(VAR$Fst,
                         as.character(VAR$seqnames), max, na.rm=TRUE))

contig.geno <- cbind(contig.geno,
                     Fis.max=tapply(VAR$Fis,
                         as.character(VAR$seqnames),
                         function(x) max(abs(x), na.rm=TRUE)))

contig.geno <- cbind(contig.geno,
                     Fis.mean=tapply(VAR$Fis,
                         as.character(VAR$seqnames), mean, na.rm=TRUE))

contig.geno <- merge(contig.geno, contig.dn.ds, by = 0)

pdf("figures/dapc.mean.dn.ds.pdf")
ggplot(contig.geno, aes(x=mean.dapc, y=contig.dn.ds)) +
    geom_point(alpha = 0.1) + geom_density2d()
dev.off()

pdf("figures/dapc.max.dn.ds.pdf")
ggplot(contig.geno, aes(x=max.dapc, y=contig.dn.ds)) +
    geom_point(alpha = 0.1) + geom_density2d()
dev.off()

pdf("figures/Fst.mean.dn.ds.pdf")
ggplot(contig.geno, aes(x=Fst.mean, y=contig.dn.ds)) +
    geom_point(alpha = 0.1) + geom_density2d()
dev.off()

pdf("figures/Fst.max.dn.ds.pdf")
ggplot(contig.geno, aes(x=Fst.max, y=contig.dn.ds)) +
    geom_point(alpha = 0.1) + geom_density2d()
dev.off()

pdf("figures/Fis.mean.dn.ds.pdf")
ggplot(contig.geno, aes(x=Fis.mean, y=contig.dn.ds)) +
    geom_point(alpha = 0.1) + geom_density2d()
dev.off()

pdf("figures/Fis.max.dn.ds.pdf")
ggplot(contig.geno, aes(x=Fis.max, y=contig.dn.ds)) +
    geom_point(alpha = 0.1) + geom_density2d()
dev.off()


double.interesting <- contig.geno[contig.geno$max.dapc > 0.00008 &
                                  contig.geno$contig.dn.ds > 0.5, "Row.names"]

MF.seperate <- TOGO.all.onto("MF", names(transcript.2.GO), 
                             double.interesting , transcript.2.GO)
GenTable(MF.seperate[[1]], MF.seperate[[2]])

BP.seperate <- TOGO.all.onto("BP", names(transcript.2.GO), 
                        double.interesting , transcript.2.GO)

GenTable(BP.seperate[[1]], BP.seperate[[2]])

CC.seperate <- TOGO.all.onto("CC", names(transcript.2.GO), 
                        double.interesting , transcript.2.GO)

GenTable(CC.seperate[[1]], CC.seperate[[2]])

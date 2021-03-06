<<vegan.pcr, echo=FALSE, cache=TRUE>>=


## this is from base
## pca.l <- lapply(contigs.l[1:6], function (x) prcomp(t(vsd[x,]), scale = FALSE, center=FALSE))
## pca.o.l <- lapply(ortho.l[1:3], function (x) prcomp(t(vsd[x,]), scale = FALSE, center=FALSE))


## samples is "sites", genes is "species"
## this here from vegan
rda.l <- lapply(contigs.l[1:6], function (x) rda(t(vsd[x,])))
rda.o.l <- lapply(ortho.l[1:3], function (x) rda(t(vsd[x,])))

rda.total <- rda(t(vsd), scale = FALSE, center=FALSE)
## netMDS <- metaMDS(t(vsd[contigs.l[[3]]),])


rda.l <- list()
rda.o.l <- list()

rda.l[[1]] <- rda(t(vsd[contigs.l[[1]], ]) ~ sex.conds, env = t(vsd[contigs.l[[1]], ]))
rda.o.l[[1]] <- rda(t(vsd[ortho.l[[1]], ]) ~ sex.conds, env = t(vsd[ortho.l[[1]], ]))

rda.l[[2]] <- rda(t(vsd[contigs.l[[2]], ]) ~ eel.conds, env = t(vsd[contigs.l[[2]], ]))
rda.o.l[[2]] <- rda(t(vsd[ortho.l[[2]], ]) ~ eel.conds, env = t(vsd[ortho.l[[2]], ]))

rda.l[[3]] <- rda(t(vsd[contigs.l[[3]], ]) ~ pop.conds, env = t(vsd[contigs.l[[3]], ]))
rda.o.l[[3]] <- rda(t(vsd[ortho.l[[3]], ]) ~ pop.conds, env = t(vsd[ortho.l[[3]], ]))


cca.l <- list()
cca.o.l <- list()

cca.l[[1]] <- cca(t(vsd[contigs.l[[1]], ]) ~ sex.conds, env = t(vsd[contigs.l[[1]], ]))
cca.o.l[[1]] <- cca(t(vsd[ortho.l[[1]], ]) ~ sex.conds, env = t(vsd[ortho.l[[1]], ]))

cca.l[[2]] <- cca(t(vsd[contigs.l[[2]], ]) ~ eel.conds, env = t(vsd[contigs.l[[2]], ]))
cca.o.l[[2]] <- cca(t(vsd[ortho.l[[2]], ]) ~ eel.conds, env = t(vsd[ortho.l[[2]], ]))

cca.l[[3]] <- cca(t(vsd[contigs.l[[3]], ]) ~ pop.conds, env = t(vsd[contigs.l[[3]], ]))
cca.o.l[[3]] <- cca(t(vsd[ortho.l[[3]], ]) ~ pop.conds, env = t(vsd[ortho.l[[3]], ]))

anosim.l <- list()
anosim.l[[1]] <- anosim(t(vsd[contigs.l[[1]],]), sex.conds,
                       permutations = 999, distance = "euclidean")
anosim.l[[2]] <- anosim(t(vsd[contigs.l[[2]],]), eel.conds,
                       permutations = 999, distance = "euclidean")
anosim.l[[3]] <- anosim(t(vsd[contigs.l[[3]],]), pop.conds,
                       permutations = 999, distance = "euclidean")

anosim.o.l <- list()
anosim.o.l[[1]] <- anosim(t(vsd[ortho.l[[1]],]), sex.conds,
                       permutations = 999, distance = "euclidean")
anosim.o.l[[2]] <- anosim(t(vsd[contigs.l[[2]],]), eel.conds,
                       permutations = 999, distance = "euclidean")
anosim.o.l[[3]] <- anosim(t(vsd[ortho.l[[3]],]), pop.conds,
                       permutations = 999, distance = "euclidean")

anova.rda.l <- lapply(rda.l, anova)
anova.rda.o.l <- lapply(rda.o.l, anova)

anova.cca.l <- lapply(cca.l, anova)
anova.cca.o.l <- lapply(cca.o.l, anova)


@

<<plot.pca, echo=FALSE, cache=TRUE>>=

## mycolors <- as.character(sex.conds:eel.conds:pop.conds)

## ## Warm colors for European population
## mycolors <- gsub("female:Aa:EU", "red", mycolors) ## female light
## mycolors <- gsub("male:Aa:EU", "darkred", mycolors) ## male dark

## mycolors <- gsub("female:Aj:EU", "yellow", mycolors)
## mycolors <- gsub("male:Aj:EU", "orange", mycolors)

## ## Cold colors for Asian population 
## mycolors <- gsub("female:Aa:TW", "green", mycolors)
## mycolors <- gsub("male:Aa:TW", "darkgreen", mycolors)

## mycolors <- gsub("female:Aj:TW", "blue", mycolors)
## mycolors <- gsub("male:Aj:TW", "darkblue", mycolors)

##  from http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual
## component 3 separates the values for population
## plot(pca.l[[3]]$x, type="n"); 
## text(pca.l[[3]]$x, rownames(pca.l[[3]]$x), cex=0.8, col=mycolors)

## library(som)
## y.som <- som(vsd[t(contigs.l[[3]]),], xdim = 5, ydim = 6, 
##              topol = "hexa", neigh = "gaussian")

## ## plot(y.som)

## library(scatterplot3d)
## library(geneplotter)

## ## smoothScatter(pca$x, col=mycolors)
## pairs(pca$x[,1:3], pch=20, col=mycolors)

## scatterplot3d(pca$x[,1:3], pch=20 , color=mycolors)

pdf("/home/ele/thesis/experimental_infection/figures/pca_sex.pdf", width=5, height=10)    
par(mfrow=c(2,1))
plot(rda.l[[1]])
plot(rda.o.l[[1]])
dev.off()

pdf("/home/ele/thesis/experimental_infection/figures/pca_eel.pdf", width=5, height=10)    
par(mfrow=c(2,1))
plot(rda.l[[2]])
plot(rda.o.l[[2]])
dev.off()

pdf("/home/ele/thesis/experimental_infection/figures/pca_pop.pdf", width=5, height=10)    
par(mfrow=c(2,1))
plot(rda.l[[3]])
plot(rda.o.l[[3]])
dev.off()

@ 

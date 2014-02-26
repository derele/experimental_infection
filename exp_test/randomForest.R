library(randomForest)

rFo.sex <- randomForest(x=t(G.e),
                        y=sex.conds,
                        ntree=1000)

rFo.eel <- randomForest(x=t(G.e),
                        y=eel.conds,
                        ntree=1000)


rFo.pop <- randomForest(x=t(G.e),
                        y=pop.conds,
                        ntree=1000)

rFo.pop.plus <- randomForest(x=cbind(t(G.e), sex.conds, eel.conds),
                             y=pop.conds,
                             ntree=1000)

rFo.eel.pop <- randomForest(x=t(G.e),
                            y=eel.conds:pop.conds,
                            ntree=1000)

rFo.all <- randomForest(x=t(G.e),
                        y=sex.conds:eel.conds:pop.conds,
                        ntree=1000)


sympa <- as.character(pop.conds:eel.conds)

sympa <- as.factor(ifelse(sympa%in%c("EU:AA", "TW:AJ"),
                          "sympa", "allopa"))

rFo.sympa <- randomForest(x=t(G.e),
                          y=sympa,
                          ntree=1000)

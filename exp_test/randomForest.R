

rFo.sex <- randomForest(x=t(G.e),
                        y=sex.conds,
                        ntree=1000)

rFo.eel <- randomForest(x=t(G.e),
                        y=eel.conds,
                        ntree=1000)


rFo.pop <- randomForest(x=t(G.e),
                        y=pop.conds,
                        ntree=1000000)

rFo.pop.plus <- randomForest(x=cbind(t(G.e), sex.conds, eel.conds),
                             y=pop.conds,
                             ntree=1000000)

rFo.eel.pop <- randomForest(x=t(G.e),
                            y=eel.conds:pop.conds,
                            ntree=10000)

rFo.all <- randomForest(x=t(G.e),
                        y=sex.conds:eel.conds:pop.conds,
                        ntree=1000)

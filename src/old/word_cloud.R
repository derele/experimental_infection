### R code from vignette source '/home/ele/thesis/experimental_infection/annot_test/word_cloud.Rnw'

###################################################
### code chunk number 1: word.cloud
###################################################


go.word.cloud <- function (terms1, terms2, lab=c("","", ""), tit=""){
  one.termtable <- terms1[, c("Term", "classic", "Significant")]
  names(one.termtable) <- c("term", "p.value.1", "freq.one")
  two.termtable <- terms2[, c("Term", "classic", "Significant")]
  names(two.termtable) <- c("term", "p.value.2", "freq.two")
  term.df<- merge(one.termtable, two.termtable, all=TRUE)
  term.df[is.na(term.df$freq.one), "freq.one"] <- 0
  term.df[is.na(term.df$freq.two), "freq.two"] <- 0
  term.df<-transform(term.df, freq.dif = freq.one-freq.two)
  term.df$p.value <- ifelse(term.df$freq.dif<0, term.df$p.value.2,
                            term.df$p.value.1)

  mi <- min(term.df$freq.dif)
  ma <- max(term.df$freq.dif)
  m <- mean(c(mi,ma))
  
  term.df <- term.df[order(term.df$freq.dif),]
  ## the occurence of the frequence difference
  occ <- rle(term.df$freq.dif)$length
  ## equal spread between min and may
  term.df$Spacing <- unlist(sapply(occ,
                                   function(x) seq(mi, ma, length=x)))
  ## patch the values where only one thing was there with mean
  term.df[term.df$Spacing==mi, "Spacing"] <- m

  ## break long names
   term.df$term <- sub("(\\w{12,}) ", "\\1\n",  term.df$term)
  term.df$p.value <- as.numeric(as.character(term.df$p.value))
  term.df$freq.one <- as.numeric(as.character(term.df$freq.one))
  
  term.cloud<-ggplot(term.df, aes(x=freq.dif, y=Spacing)) +
    geom_text(aes(size=freq.one, label=term, colour=p.value)) +
      scale_size(to=c(3,11),
                 name="GO-annotation Frequency") +
                   scale_colour_gradient("p-value for overrepresentation",
                                         breaks=c(0.001, 0.01, 0.02, 0.03, 0.04),
                                         low="darkred",
                                         high="darkblue"
                                         )
  ## aesthetic corrections
  term.cloud <- term.cloud +
    scale_x_continuous(breaks=c(mi, 0, ma), labels=rev(lab))

  ## aesthetic corrections
  term.cloud <- term.cloud + 
    scale_y_continuous(breaks=c(0),
                       labels=c("")) +
                         xlab("") +
                           ylab("")+
                             theme_bw()+
                               opts(panel.grid.major=theme_blank(),
                                    panel.grid.minor=theme_blank(),
                                    title=tit)
 return(term.cloud)
}


go.eel.vs.pop.mf <- go.word.cloud(GO.over.l[[3]],
                                  GO.over.l[[2]],
                                  c("Enriched more in population differences",
                                    "\nEnriched equally",
                                    "Enriched more in host-species differences"),
                                  "Frequency of enriched GO-terms for nematode-population vs. host-species")

ggsave("/home/ele/thesis/experimental_infection/figures/word_cloud_eel_vs_pop_mf.pdf", go.eel.vs.pop.mf, width=22, height=8)




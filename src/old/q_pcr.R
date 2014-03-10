library(ggplot2)
library(reshape)

###########################################
### Cox ###################################

COX <- read.delim("/home/ele/thesis/experimental_infection/q_pcr/COX_Rad.CSV", sep="," )

mean.techrep <- function (x) tapply(x$Ct.SYBR,
                                    list(x$Target.SYBR, x$Name),
                                    mean, na.rm=TRUE)

COX <- mean.techrep(COX)


exp.delta.ct <- function(x, norm) {
  d.ct <- apply(x, 1, function (y)
                y - x[grepl(norm, rownames(x)),])
  e.d.ct <- as.data.frame(2^-(d.ct))
  e.d.ct[, !grepl(norm, names(e.d.ct))]
}
  
COX.dct <- exp.delta.ct(COX, "60S")


COX.dct$eel <- ifelse(grepl("R|UW", rownames(COX.dct)), "Euro", "Asia")
COX.dct$worm <- rownames(COX.dct)

COX <- melt(COX.dct)

ggplot(COX, aes(eel, value, label=worm)) + geom_text() +
  facet_wrap(~variable, scales="free") 


##########################################
### Sex ##################################


## SEX <- read.delim("/home/ele/thesis/experimental_infection/q_pcr/Sex1.CSV", sep=",")

SEX1 <- read.delim("/home/ele/thesis/experimental_infection/q_pcr/Sex1.CSV", sep=",")
SEX1 <- mean.techrep(SEX1)
SEX1.dct <- exp.delta.ct(SEX1, "60S")

SEX2 <- read.delim("/home/ele/thesis/experimental_infection/q_pcr/Sex2.CSV", sep=",")
SEX2 <- mean.techrep(SEX2)
SEX2.dct <- exp.delta.ct(SEX2, "60S")
SEX2.dct$"7271qs" <- NA


SEX3 <- read.delim("/home/ele/thesis/experimental_infection/q_pcr/Sex3.CSV", sep=",")
SEX3 <- mean.techrep(SEX3)
SEX3.dct <- exp.delta.ct(SEX3, "60S")

SEX2 <- SEX2[, names(SEX1)]
SEX3 <- SEX3[, names(SEX1)]

SEX <- rbind(SEX1.dct,
             SEX2.dct,
             SEX3.dct)

SEX$sex <- ifelse(grepl("M", rownames(SEX)), "male", "female")
SEX$worm <- rownames(SEX)

SEX <- melt(SEX, direction="long")

ggplot(SEX, aes(sex, value, label=worm)) + geom_text() +
  facet_wrap(~variable, scales="free") 

## looks very bad from Contig76 and Contig86, which should be high in
## females, low in males. It is only one female however. Look at the
## raw data and decide if to exclude!!

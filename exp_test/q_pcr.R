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
  as.data.frame(2^-(d.ct))}
  
COX.dct <- exp.delta.ct(COX, "60S")


COX.dct$eel <- ifelse(grepl("R|UW", rownames(COX.dct)), "Euro", "Asia")


ggplot(COX.dct, aes(eel, COX2.1)) + geom_boxplot()

##########################################
### Sex ##################################


## SEX <- read.delim("/home/ele/thesis/experimental_infection/q_pcr/Sex1.CSV", sep=",")

SEX2 <- read.delim("/home/ele/thesis/experimental_infection/q_pcr/Sex2.CSV", sep=",")
SEX2 <- mean.techrep(SEX2)
SEX2.dct <- exp.delta.ct(SEX2, "60S")


SEX1 <- read.delim("/home/ele/thesis/experimental_infection/q_pcr/Sex1.CSV", sep=",")
SEX1 <- mean.techrep(SEX1)
SEX1.dct <- exp.delta.ct(SEX1, "60S")


SEX <- rbind(SEX1.dct[, c(1, 3, 5, 6)],
             SEX2.dct[, c(1, 3, 4, 5)])

SEX$sex <- ifelse(grepl("M", rownames(SEX)), "male", "female")
SEX <- melt(SEX, direction="long")

ggplot(SEX, aes(sex, value)) + geom_jitter() + facet_wrap(~variable, scales="free")

## looks very bad from Contig76 and Contig86, which should be high in
## females, low in males. It is only one female however. Look at the
## raw data and decide if to exclude!!

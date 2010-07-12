#Emanuel Heitilinger Basic analysis of experimental infections 
EA <- read.delim("AAL_SCHL_EH.csv", header=TRUE, sep=",", as.is=TRUE)

facs <- c("A.No.", "A.spec", "L.population", "Dpi","tank","L2")
measure <- c("length", "weight")
counts <- c("l3", "l4", "adult.m", "adult.f","enc")
ofacs <- c("egg","Pseudo.0.3")

EA[,facs] <- lapply(EA[,facs], as.factor)
EA[,measure] <- lapply(EA[,measure], as.numeric)
EA[,counts] <- lapply(EA[,counts], as.numeric)
EA[,ofacs] <- lapply(EA[,ofacs],function (x) as.ordered(as.numeric(as.character(x))))

E <- EA[, c(facs, measure, counts, ofacs)]
E$intensity <- apply(data.frame(E$adult.m, E$adult.f), 1, sum)

c(by(E, E$A.spec:E$L.population, function (x) nrow(x[x$adult.m>0,])))
c(by(E, E$A.spec:E$L.population, function (x) nrow(x[x$adult.f>0,])))
c(by(E, E$A.spec:E$L.population, function (x) sum(x$adult.m)))
c(by(E, E$A.spec:E$L.population, function (x) sum(x$adult.f)))

# use later 
dropfa <- function (x) factor(x, levels=unique(x))

AA.T.F <- subset(E, E$L.population=="T" & E$A.spec=="AA" & E$adult.f>0)
AA.T.M <- subset(E, E$L.population=="T" & E$A.spec=="AA" & E$adult.m>0)
AA.R.M <- subset(E, E$L.population=="R" & E$A.spec=="AA" & E$adult.m>0)
AA.R.F <- subset(E, E$L.population=="R" & E$A.spec=="AA" & E$adult.f>0)
AJ.R.F <- subset(E, E$L.population=="R" & E$A.spec=="AJ" & E$adult.f>0)
AJ.R.M <- subset(E, E$L.population=="R" & E$A.spec=="AJ" & E$adult.m>0)
AJ.T.M <- subset(E, E$L.population=="T" & E$A.spec=="AJ" & E$adult.m>0)
AJ.T.F <- subset(E, E$L.population=="T" & E$A.spec=="AJ" & E$adult.f>0)


# The three biological replicates will be chosen to contain one worm
# from a low infection (intensity==1),
# one from a medium infection (intesity==round(mean(intensity))) and
# one from a high infection (intesity==max(intensity))


#################### CHOOSE WORMS for PREP ###################### 
#prep <- as.data.frame()
#names(prep) <- names(E)
set.seed(123) # set random seed
# Female Anguilla anguilla, Taiwanese


prep <- as.data.frame(
              AA.T.F[sample(1:nrow(AA.T.F),1),])
t(prep[nrow(prep),])
# so the first worm RNA-prepped is AA/T-20 it has intensity 1
# WORKED

#say what sex was prepped:
#prep <- cbind(sex="F", prep)

ME <- round(mean(AA.T.F$intensity))
prep <- rbind(prep, 
              AA.T.F[AA.T.F$intensity==ME,]
              [sample(1:nrow(AA.T.F[AA.T.F$intensity==ME,]),1),])
t(prep[nrow(prep),])
# AA/T-32
# FAILED low concentration
prep <- rbind(prep,
              AA.T.F[AA.T.F$intensity==ME,]
              [sample(1:nrow(AA.T.F[AA.T.F$intensity==ME,]),1),])
t(prep[nrow(prep),])
# AA/T-45 (one female, so try another one to be sure)
prep <- rbind(prep,
              AA.T.F[AA.T.F$intensity==ME,]
              [sample(1:nrow(AA.T.F[AA.T.F$intensity==ME,]),1),])
t(prep[nrow(prep),])
#AA/T-46 3 females
prep <- rbind(prep,
              AA.T.F[AA.T.F$intensity==ME,]
              [sample(1:nrow(AA.T.F[AA.T.F$intensity==ME,]),1),])
t(prep[nrow(prep),])
#AA/T-1 two females

MAX <- max(AA.T.F$intensity)
prep <- rbind(prep,
              AA.T.F[AA.T.F$intensity==MAX,]
              [sample(1:nrow(AA.T.F[AA.T.F$intensity==MAX,]),1),])
t(prep[nrow(prep),])
# prep from AA/T-12
# WORKED

#prep$sex.prepped <- "female"


# Male Anguilla anguilla, Taiwanese
# tank 9 and 17 preferred, becaus we have
# no female samples from those

prep <- rbind(prep, AA.T.M[(AA.T.M$tank==9|AA.T.M$tank==17) &
              AA.T.M$intensity==1,]
              [sample(1:nrow(AA.T.M[(AA.T.M$tank==9|AA.T.M$tank==17)
              &AA.T.M$intensity==1,]),1),])
t(prep[nrow(prep),])
# AA/T-36 but this is the only one with intensity==1
# one with intensity 2 to make sure to get a low intensity sample
prep <- rbind(prep, AA.T.M[(AA.T.M$tank==9|AA.T.M$tank==17) &
              AA.T.M$intensity==2,]
              [sample(1:nrow(AA.T.M[(AA.T.M$tank==9|AA.T.M$tank==17)
              &AA.T.M$intensity==2,]),1),])
t(prep[nrow(prep),])
# AA/T40 one male
prep <- rbind(prep, AA.T.M[(AA.T.M$tank==9|AA.T.M$tank==17) &
              AA.T.M$intensity==2,]
              [sample(1:nrow(AA.T.M[(AA.T.M$tank==9|AA.T.M$tank==17)
              &AA.T.M$intensity==2,]),1),])
t(prep[nrow(prep),])
# two males AA/T38

ME <- round(mean(AA.T.M$intensity))


prep <- rbind(prep, AA.T.M[(AA.T.M$tank==9|AA.T.M$tank==17) &
              AA.T.M$intensity==ME,]
              [sample(1:nrow(AA.T.M[(AA.T.M$tank==9|AA.T.M$tank==17)
              &AA.T.M$intensity==ME,]),1),])
t(prep[nrow(prep),])
# AA/T-43 4 males now enogh of tank 17

prep <- rbind(prep, AA.T.M[(AA.T.M$tank==9) &
                           AA.T.M$intensity==ME+1,]
              [sample(1:nrow(AA.T.M[(AA.T.M$tank==9)
                                    &AA.T.M$intensity==ME+1,]),1),])
t(prep[nrow(prep),])
# AA/T-24 4 males

#The maimal infected is already prepped for female
MAX <- 11 # this is the second highest

prep <- rbind(prep, AA.T.M[AA.T.M$intensity==MAX,]
              [sample(1:nrow(AA.T.M[AA.T.M$intensity==MAX,]),1),])
t(prep[nrow(prep),])

# test, that we have no eel twice
length(prep$A.No.)==length(unique(prep$A.No.))


 ####
# Female Anguilla japonica, Taiwanese

# intensity is 1 and a female is found
prep <- rbind(prep, AJ.T.F[AJ.T.F$adult.f==1 & AJ.T.F$intensity==1,]
              [sample(1:nrow(AJ.T.F[AJ.T.F$adult.f==1 & AJ.T.F$intensity==1,]),1),])
t(prep[nrow(prep),])
# AJ/T7 two other low intensity ones
prep <- rbind(prep, AJ.T.F[AJ.T.F$adult.f>0 & AJ.T.F$intensity<3,]
              [sample(1:nrow(AJ.T.F[AJ.T.F$adult.f>0 & AJ.T.F$intensity<3,]),1),])
t(prep[nrow(prep),])
# AJ/T-3

prep <- rbind(prep, AJ.T.F[AJ.T.F$adult.f>0 & AJ.T.F$intensity<3,]
              [sample(1:nrow(AJ.T.F[AJ.T.F$adult.f>0 & AJ.T.F$intensity<3,]),1),])
t(prep[nrow(prep),])
# AJ/T-26

ME <- round(mean(AJ.T.F$intensity))

prep <- rbind(prep, AJ.T.F[AJ.T.F$intensity>ME-1 & AJ.T.F$intensity<ME+1,]
              [sample(1:nrow(AJ.T.F[AJ.T.F$intensity>ME-1 & AJ.T.F$intensity<ME+1,]),1),])
t(prep[nrow(prep),])
# AJ/T-8

MAX <- max(AJ.T.F$intensity)

# there are two eels with this value so the random choosing will
# work also between males and females.
prep <- rbind(prep, AJ.T.F[AJ.T.F$intensity==MAX,]
              [sample(1:nrow(AJ.T.F[AJ.T.F$intensity==MAX,]),1),])
t(prep[nrow(prep),])
# AJ/T19

# Repeating one low intensitysample, because first was rather low
AJ.T.F[order(AJ.T.F$intensity),]
prep <- rbind(prep, AJ.T.F[AJ.T.F$A.No.=="AJ7T-4"|
                           AJ.T.F$A.No.=="AJ/T-5",])


##############
# Non randomly sampled low number of worms:
# AJ/R put in here!
##############

# Male Anguilla japonica, Taiwanese
# Pick by looking at the sorted values

#AJ/T12

AJ.T.M[order(AJ.T.M$intensity),]

#low intesity samples
#AJ/T-12 AJ/T-4 AJ/T-24 AJ/T-6
prep <- rbind(prep, AJ.T.M[AJ.T.M$A.No.=="AJ/T-4"|
                           AJ.T.M$A.No.=="AJ/T-24"|
                           AJ.T.M$A.No.=="AJ/T-6",])

#medium intensity samples
#AJ/T10 AJ/T23
prep <- rbind(prep, AJ.T.M[AJ.T.M$A.No.=="AJ/T-10",])

#High intensity
#AJ/T-25
prep <- rbind(prep, AJ.T.M[AJ.T.M$A.No.=="AJ/T-25",])

# For all low and medium intensities preps did not work!
# Picking anoterh 6 samples
prep <- rbind(prep, AJ.T.M[AJ.T.M$A.No.=="AJ/-T-30"|
                           AJ.T.M$A.No.=="AJ/T-20"|
                           AJ.T.M$A.No.=="AJ/T-27",])

###############################################
# samples from European eels Japanese larvae females
# AA/R-F
AA.R.F[order(AA.R.F$intensity),]

prep <- rbind(prep, AA.R.F[AA.R.F$A.No.=="AA/R-29"|
                           AA.R.F$A.No.=="AA/R-18"|
                           AA.R.F$A.No.=="AA/R-2"|
                           AA.R.F$A.No.=="AA/R-28"|
                           AA.R.F$A.No.=="AA/R-8",])


AA.R.M[order(AA.R.M$intensity),]

prep <- rbind(prep, AA.R.M[AA.R.M$A.No.=="AA/R-3"|
                           AA.R.M$A.No.=="AA/R-25"|
                           AA.R.M$A.No.=="AA/R-4"|
                           AA.R.M$A.No.=="AA/R-16"|
                           AA.R.M$A.No.=="AA/R-11",])




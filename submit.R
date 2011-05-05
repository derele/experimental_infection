# This file is a indirect result of eeldata.R. In the latter I
# selected female worms randomly and male worms accordingly not at
# random (because of size and preparation issues the biggest worms
# were used

# The following lines just select what is handed in to the GenePool
# for tag-sequencing



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

## AA/T F

SU <- cbind(E[E$A.No=="AA/T-20",], prepsex="female", prepnum=1,
            prepconc=5.6)

SU <- rbind(SU, cbind(E[E$A.No=="AA/T-12",], prepsex="female",
                     prepnum=1, prepconc=6.8))

SU <- rbind(SU, cbind(E[E$A.No=="AA/T-45",], prepsex="female",
                     prepnum=1, prepconc=8))

## AA/T M

SU <- rbind(SU, cbind(E[E$A.No=="AA/T-24",], prepsex="male", prepnum=3,
            prepconc=4.8))

SU <- rbind(SU, cbind(E[E$A.No=="AA/T-42",], prepsex="male", prepnum=1,
            prepconc=5.6))

SU <- rbind(SU, cbind(E[E$A.No=="AA/T-3",], prepsex="male", prepnum=4,
            prepconc=4.88))

## AA/R F

SU <- rbind(SU, cbind(E[E$A.No=="AA/R-18",], prepsex="female",
                     prepnum=1, prepconc=4.8))

SU <- rbind(SU, cbind(E[E$A.No=="AA/R-28",], prepsex="female",
                     prepnum=1, prepconc=5.2))

SU <- rbind(SU, cbind(E[E$A.No=="AA/R-8",], prepsex="female",
                     prepnum=1, prepconc=5.2))

## AA/R M

SU <- rbind(SU, cbind(E[E$A.No=="AA/R-16",], prepsex="male",
                     prepnum=4, prepconc=5.2))

SU <- rbind(SU, cbind(E[E$A.No=="AA/R-11",], prepsex="male",
                     prepnum=14, prepconc=6.4))

SU <- rbind(SU, cbind(E[E$A.No=="AA/R-2",], prepsex="male",
                     prepnum=4, prepconc=6.6))

## AJ/T F

SU <- rbind(SU, cbind(E[E$A.No=="AJ/T-8",], prepsex="female",
                     prepnum=1, prepconc=5.91))

SU <- rbind(SU, cbind(E[E$A.No=="AJ/T-5",], prepsex="female",
                     prepnum=1, prepconc=4.8))

### NOT SURE ABOUT THIS CONC SOO LOW

SU <- rbind(SU, cbind(E[E$A.No=="AJ/T-4",], prepsex="female",
                     prepnum=1, prepconc=2.8))

## AJ/T M

SU <- rbind(SU, cbind(E[E$A.No=="AJ/T-25",], prepsex="male",
                     prepnum=5, prepconc=4.05))

SU <- rbind(SU, cbind(E[E$A.No=="AJ/T-19",], prepsex="male",
                     prepnum=7, prepconc=3.5))

SU <- rbind(SU, cbind(E[E$A.No=="AJ/T-20",], prepsex="male",
                     prepnum=8, prepconc=3.8))
## AJ/R F

SU <- rbind(SU, cbind(E[E$A.No=="AJ/R-1",], prepsex="female",
                     prepnum=1, prepconc=5.92))

SU <- rbind(SU, cbind(E[E$A.No=="AJ/R-3",], prepsex="female",
                     prepnum=1, prepconc=6.9))

SU <- rbind(SU, cbind(E[E$A.No=="AJ/R-5",], prepsex="female",
                     prepnum=1, prepconc=4.04))
## AJ/R M

SU <- rbind(SU, cbind(E[E$A.No=="AJ/R-1",], prepsex="male",
                     prepnum=1, prepconc=2.5))

SU <- rbind(SU, cbind(E[E$A.No=="AJ/R-3",], prepsex="male",
                     prepnum=2, prepconc=2.6))

SU <- rbind(SU, cbind(E[E$A.No=="AJ/R-5",], prepsex="male",
                     prepnum=1, prepconc=2.23))

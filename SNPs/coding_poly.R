## ## BASE ONTOLOGY
## #
## # Code from Mark Blaxter, modified by John Davey,
## # Translated from Perl to R by Emanuel Heitlinger:

## # A phase 1 any change is nonsynonymous
## # B phase 2 any change is nonsynonymous
## # C phase 3 any change is nonsynonymous
## # D phase 1 change to CT is nonsynonymous
## # E phase 2 change to CT is nonsynonymous
## # F phase 3 change to CT is nonsynonymous
## # G phase 1 change to AG is nonsynonymous
## # H phase 2 change to AG is nonsynonymous
## # I phase 3 change to AG is nonsense
## # K phase 1 change to GT is nonsynonymous
## # L phase 2 change to A is nonsense, to anything else is nonsynonymous
## # J phase 3 change to G is nonsynonymous
## # M phase 3 change to G is nonsense, to A is nonsynonymous
## # N phase 3 any change synonymous
## # O phase 1 change to T nonsense, others nonsynonymous
## # P phase 3 change to AG is nonsynonymous
## # Q phase 1 change to T nonsense, to G nonsynonymous
## # R phase 2 change to AG nonsense, others nonsynonymous
## # S phase 3 change to A nonsense, others nonsynonymous
## # T phase 3 change to A nonsense, G nonsynonymous

## # W all changes are unknown # EH added 08/23/2011

## #        a           g           c           t
## #
## # a     aaa K OBF   aga R QBF   aca T ABN   ata I ABJ
## #       aag K OBF   agg R KBF   acg T ABN   atg M ABC
## #       aac N ABP   agc S ABP   acc T ABN   atc I ABJ
## #       aat N ABP   agt S ABP   act T ABN   att I ABJ
## #
## # g     gaa E OBF   gga G OBN   gca A ABN   gta V ABN
## #       gag E OBF   ggg G ABN   gcg A ABN   gtg V ABN
## #       gac D ABP   ggc G ABN   gcc A ABN   gtc V ABN
## #       gat D ABP   ggt G ABN   gct A ABN   gtt V ABN
## #
## # c     caa Q OBF   cga R QBN   cca P ABN   cta L GBN
## #       cag Q OBF   cgg R KBN   ccg P ABN   ctg L GBN
## #       cac H ABP   cgc R ABN   ccc P ABN   ctc L ABN
## #       cat H ABP   cgt R ABN   cct P ABN   ctt L ABN
## #
## # t     taa * AEF   tga * AEC   tca S ARN   tta L GRF
## #       tag * ABF   tgg W ALS   tcg S ALN   ttg L GLF
## #       tac Y ABI   tgc C ABT   tcc S ABN   ttc F ABP
## #       tat Y ABI   tgt C ABT   tct S ABN   ttt F ABP





    ## # Set up Base Ontology vector
base.ontology.encode.string = c(
  "aaa" = "OBF",
  "aag" = "OBF",
  "aac" = "ABP",
  "aat" = "ABP",
  "aga" = "QBF",
  "agg" = "KBF",
  "agc" = "ABP",
  "agt" = "ABP",
  "aca" = "ABN",
  "acg" = "ABN",
  "acc" = "ABN",
  "act" = "ABN",
  "ata" = "ABJ",
  "atg" = "ABC",
  "atc" = "ABJ",
  "att" = "ABJ",
  ##
  "gaa" = "OBF",
  "gag" = "OBF",
  "gac" = "ABP",
  "gat" = "ABP",
  "gga" = "OBN",
  "ggg" = "ABN",
  "ggc" = "ABN",
  "ggt" = "ABN",
  "gca" = "ABN",
  "gcg" = "ABN",
  "gcc" = "ABN",
  "gct" = "ABN",
  "gta" = "ABN",
  "gtg" = "ABN",
  "gtc" = "ABN",
  "gtt" = "ABN",
  ##
  "caa" = "OBF",
  "cag" = "OBF",
  "cac" = "ABP",
  "cat" = "ABP",
  "cga" = "QBN",
  "cgg" = "KBN",
  "cgc" = "ABN",
  "cgt" = "ABN",
  "cca" = "ABN",
  "ccg" = "ABN",
  "ccc" = "ABN",
  "cct" = "ABN",
  "cta" = "GBN",
  "ctg" = "GBN",
  "ctc" = "ABN",
  "ctt" = "ABN",
  ##
  "taa" = "AEF",
  "tag" = "ABF",
  "tac" = "ABI",
  "tat" = "ABI",
  "tga" = "AEC",
  "tgg" = "ALS",
  "tgc" = "ABT",
  "tgt" = "ABT",
  "tca" = "ARN",
  "tcg" = "ALN",
  "tcc" = "ABN",
  "tct" = "ABN",
  "tta" = "GRF",
  "ttg" = "GLF",
  "ttc" = "ABP",
  "ttt" = "ABP"
  );

base.ontology.encode <- function(x){  
  ## # wrap in a control function for iupac and other "bad"
  ## # bases
  if (nchar(gsub("[bdefhijklmnopqrsuvwxyz]", "", x, ignore.case=TRUE)) != 3){
    return (paste(rep("W", times = nchar(x)), collapse=""))
  }
  else {
    return(base.ontology.encode.string[x])
  }
}
  
base.ontology.decode = list(
  "A" = c("A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "B" = c("A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "C" = c( "A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "D" = c( "A" = "Synonymous", "C" = "Nonsynonymous",
    "G" = "Synonymous", "T" = "Nonsynonymous" ),
  "E" = c( "A" = "Synonymous", "C" = "Nonsynonymous",
    "G" = "Synonymous", "T" = "Nonsynonymous" ),
  "F" = c( "A" = "Synonymous", "C" = "Nonsynonymous",
    "G" = "Synonymous", "T" = "Nonsynonymous" ),
  "G" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "H" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "I" = c( "A" = "Nonsense",  "C" = "Synonymous",
    "G" = "Nonsense",  "T" = "Synonymous" ),
  "J" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "K" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "L" = c( "A" = "Nonsense",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "M" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsense",  "T" = "Synonymous" ),
  "N" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Synonymous", "T" = "Synonymous" ),
  "O" = c( "A" = "Nonsynonymous",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsense" ),
  "P" = c( "A" = "Nonsynonymous",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "Q" = c( "A" = "Synonymous", "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsense" ),
  "R" = c( "A" = "Nonsense",  "C" = "Nonsynonymous",
    "G" = "Nonsense",  "T" = "Nonsynonymous" ),
  "S" = c( "A" = "Nonsense",  "C" = "Nonsynonymous",
    "G" = "Nonsynonymous",  "T" = "Nonsynonymous" ),
  "T" = c( "A" = "Nonsense",  "C" = "Synonymous",
    "G" = "Nonsynonymous",  "T" = "Synonymous" ),
  "X" = c( "A" = "Nonsense",  "C" = "Nonsense",
    "G" = "Nonsense",  "T" = "Nonsense" ),
  "W" = c("A" = NA,  "C" = NA,
    "G" = NA, "T" = NA , "R" = NA,
    "Y" = NA, "S" = NA, "W" = NA,
    "K" = NA, "M" = NA, "B" = NA,
    "D" = NA, "H" = NA, "V" = NA,
    "N" = NA, "X" = NA),
  "Y" = c("A" = "low coverage",  "C" = "low coverage",
    "G" = "low coverage", "T" = "low coverage" , "R" = "low coverage",
    "Y" = "low coverage", "S" = "low coverage", "W" = "low coverage",
    "K" = "low coverage", "M" = "low coverage", "B" = "low coverage",
    "D" = "low coverage", "H" = "low coverage", "V" = "low coverage",
      "N" = "low coverage", "X" = "low coverage"),
  "Z" = c("A" = "outside ORF",  "C" = "outside ORF",
    "G" = "outside ORF", "T" = "outside ORF" , "R" = "outside ORF",
    "Y" = "outside ORF", "S" = "outside ORF", "W" = "outside ORF",
    "K" = "outside ORF", "M" = "outside ORF", "B" = "outside ORF",
    "D" = "outside ORF", "H" = "outside ORF", "V" = "outside ORF",
    "N" = "outside ORF", "X" = "outside ORF")
  );


## read the gff of protein annotation
library(rtracklayer)

transcripts <- read.ESeq("/data/RNAseq/protein_prediction/Trinity.fasta")

names(transcripts) <- gsub("(^comp\\d+_c\\d+_seq\\d+) +.*",
                           "\\1",
                           names(transcripts))

transcripts.1 <- transcripts[gsub(".*_(seq\\d+)", "\\1",
                                  names(transcripts))%in%"seq1"]


all.gff <- import.gff("/data/RNAseq/protein_prediction/best_candidates.gff3",
                       asRangedData=FALSE)

cds.gff <- subset(all.gff, all.gff$type=="CDS")
cds.1.gff <- subset(cds.gff, as.character(cds.gff@seqnames)%in%names(transcripts.1))
cds.1.gff.df <- as.data.frame(cds.1.gff)
cds.1.gff.df <- merge(cds.1.gff.df, transcripts.1,
                      by.x="seqnames",
                      by.y=0)
names(cds.1.gff.df)[names(cds.1.gff.df)%in%"y"] <- "transcript"

cds <- apply(cds.1.gff.df, 1, function(x) {
  cds <- substr(x["transcript"],
                x["start"],
                x["end"])
  if(x["strand"]%in%"+"){
    return(cds)
  }
  if(x["strand"]%in%"-"){
    revcom(cds)
  }
})

names(cds) <- cds.1.gff.df$seqnames
cds.1.gff.df$cds <- cds


get.ontology <- function (cds.seq){
  codons <- substring(cds.seq,
                      seq(1, nchar(cds.seq), by=3),
                      seq(3,nchar(cds.seq), by=3))
  ontology <- lapply(tolower(codons),  base.ontology.encode)
  return(paste(ontology, collapse=""))
}

base.ontology.cds <- sapply(cds, get.ontology)
names(base.ontology.cds) <- names(cds)

cds.1.gff.df$base.ontology.cds <- base.ontology.cds


foo <- apply(cds.1.gff.df, 1, function (x) {
  if(x["strand"]%in%"-"){
    ont <- strReverse(x["base.ontology.cds"])
  }
   if(x["strand"]%in%"+"){
     ont <- x["base.ontology.cds"]
  }
  start <- as.numeric(as.character(x["start"]))
  end <- as.numeric(as.character(x["end"]))
  before <- ifelse( start > 1, start-1, 0)
  after <- nchar(as.character(x["transcript"])) - end + 1
  paste(paste(rep("Z", times = before ), collapse=""),
        ont,
        paste(rep("Z", times = after ), collapse=""),
        sep = "")
})


bad <- (nchar(as.character(cds.1.gff.df$transcript)) -
        as.numeric(as.character(cds.1.gff.df$end)) +1)<0

num <- nchar(as.character(cds.1.gff.df$transcript)) - as.numeric(as.character(cds.1.gff.df$end))

summary.factor(bad&bad2)

source("functions.R")

if(!exists("VCF")){
    source("./SNPs/genotyping.R")
}

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

base.ontology.encode <- function(x){  
    ## Set up Base Ontology vector
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

transcripts <- read.ESeq("/data/A_crassus/RNAseq/Trinity/Trinity.fasta")
names(transcripts) <- gsub("(^comp\\d+_c\\d+_seq\\d+) +.*",
                           "\\1",
                           names(transcripts))

transcripts.1 <- transcripts[gsub(".*_(seq\\d+)", "\\1",
                                  names(transcripts))%in%"seq1"]


all.gff <- import.gff("/data/A_crassus/RNAseq/protein_prediction/best_candidates.gff3",
                       asRangedData=FALSE)

cds.gff <- subset(all.gff, all.gff$type=="CDS")
cds.1.gff <- subset(cds.gff, as.character(cds.gff@seqnames)%in%
                    names(transcripts.1))
cds.1.gff.df <- as.data.frame(cds.1.gff)
cds.1.gff.df <- merge(cds.1.gff.df, transcripts.1,
                      by.x="seqnames",
                      by.y=0)

names(cds.1.gff.df)[names(cds.1.gff.df)%in%"y"] <- "transcript"
cds.1.gff.df$transcript <- as.character(cds.1.gff.df$transcript)
cds.1.gff.df$start <- as.numeric(cds.1.gff.df$start)
cds.1.gff.df$end <- as.numeric(cds.1.gff.df$end)

get.ontology <- function(transcript, start, end, strand) {
    utr1 <- substr(transcript, 1,
                    as.numeric(start) - 1)
    coding <- tolower(substr(transcript, as.numeric(start),
                             as.numeric(end)))
    if(nchar(coding)%%3!=0){
        warning("coding region has not multiple of 3 lenght:\n",
                transcript, "\tstart:", start,
                "\tend:", end, "\tstrand:", strand)
    }
    utr2 <- substr(transcript,
                    as.numeric(end)+1,
                    nchar(transcript))
    if(strand%in%"-"){
        coding <- revcom(coding) # revcom thins on the minus strand
    }
    ## split the cds by 3 and get the base ontology for each codon
    codons <- substring(coding,
                        seq(1, nchar(coding), by=3),
                        seq(3,nchar(coding), by=3))
    ont <- lapply(codons, base.ontology.encode)    
    ontology <- paste(ont, collapse="")
    utr1 <- gsub("\\w", "Z", utr1)
    utr2 <- gsub("\\w", "Z", utr2)
    if(strand%in%"+"){
        return(paste(utr1, ontology, utr2, sep=""))
    }
    ### Need to strReverse the wrong way round cds
    if(strand%in%"-"){
        return(paste(utr1, strReverse(ontology), utr2, sep=""))
    }
}

cds.1.gff.df$base.ontology <-
    as.character(apply(cds.1.gff.df, 1, function (x) {
        get.ontology(x["transcript"], x["start"],
                     x["end"], x["strand"])}))

## all ontologies have the right length
summary.factor(nchar(cds.1.gff.df$transcript)==
               nchar(cds.1.gff.df$base.ontology))

VCF$SNP <- rownames(VCF)
VAR <- merge(cds.1.gff.df, VCF, by.x = "seqnames", by.y = "V1")

get.effect <- function (ontology, base, to, strand){
    if(strand%in%"-"){
        to <- revcom(to)
    }
    code <- unlist(strsplit(ontology, ""))[[base]]
    base.ontology.decode[[code]][to]
}

VAR$effect <- apply(VAR, 1, function(x){
  get.effect(x["base.ontology"], as.numeric(x["V2"]), x["V5"],
             as.character(x["strand"]))
})

## somehow the minus strand has much higher N/S
tapply(VAR$effect, as.character(VAR$strand), table)

## but also SNPs seperating populations have are no more likely to be Nonsyms 
tapply(VAR$effect, as.character(VAR$SNP)%in%rownames(GT.haplo.diff), table)

get.sites <- function (ontology, transcript){
  s.sites <- sapply(1:length(ontology), function (i) {
     split.ont <- unlist(strsplit(ontology[[i]], ""))
     split.cod <- unlist(strsplit(transcript[[i]], ""))
     decoded <- lapply(split.ont, function (a){
       base.ontology.decode[[a]]})
     reduced.decoded <- lapply(1:length(decoded), function (x) {
       subset(decoded[[x]], names(decoded[[x]])!=toupper(split.cod[[x]]))})
     s <- lapply(reduced.decoded, function (w) {
       length(w[grepl("Syn", w)])/length(w[grepl("Syn|Non", w)])})
     n <- lapply(reduced.decoded, function (w) {
       length(w[grepl("Non", w)])/length(w[grepl("Syn|Non", w)])})
     nsyn.sites <- sum(unlist(n), na.rm=TRUE)
     syn.sites <- sum(unlist(s), na.rm=TRUE)
     cbind(nsyn.sites, syn.sites)
   })
  data.frame(t(s.sites))
}

sites <- get.sites(VAR[!duplicated(VAR$seqnames), "base.ontology"],
                   VAR[!duplicated(VAR$seqnames), "transcript"])

rownames(sites) <-  VAR$seqnames[!duplicated(VAR$seqnames)]
names(sites) <- c("nsyn.sites", "syn.sites")

VAR <- merge(VAR, sites, by.x = "seqnames", by.y = 0)

get.dn.ds <- function(VARobj){
    uni <- !duplicated(as.character(VARobj[,"seqnames" ]))
    (nrow(VARobj[grepl("Non*", VARobj$effect),])/
     sum(as.numeric(VARobj[uni, "nsyn.sites"]), na.rm=T))/
     (nrow(VARobj[VARobj$effect=="Synonymous",])/
      sum(as.numeric(VARobj[uni, "syn.sites"]), na.rm=T))
}

dn.ds.overall <- get.dn.ds(VAR)

## there is something severely wrong with the minus strand!!!
by(VAR, as.character(VAR$strand), get.dn.ds)

contig.dn.ds <- by(VAR, as.character(VAR$seqnames), get.dn.ds)

contig.dn.ds[is.infinite(contig.dn.ds)] <- 4
contig.dn.ds[is.na(contig.dn.ds)] <- 0

pos.selected <- (names(contig.dn.ds[contig.dn.ds>0.5 &
                                    contig.dn.ds<5]))

## SNPs unique to a population are producing higher dn/ds
by(VAR, as.factor(as.character(VAR$SNP)%in%rownames(GT.haplo.diff)):
   as.factor(as.character(VAR$strand)), get.dn.ds)

sapply(1:11, function (i) {
    get.dn.ds(VAR[VAR$SNP%in%GT.haplo.diff[[i]], ])
})

sapply(1:11, function (i) {
    get.dn.ds(VAR[VAR$SNP%in%GT.haplo.diff[[i]], ])
})


PCA.loadings <- pca1$loadings
rownames(PCA.loadings) <- locNames(GT.ade)

VAR$SNP <- gsub("\\.", "_" , VAR$SNP)

VAR <- merge(VAR, PCA.loadings, by.x = "SNP", by.y = 0)

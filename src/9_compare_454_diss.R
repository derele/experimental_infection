library(reshape)
library(limma)

if(!exists("read.ESeq")){
    source("src/1_functions.R")
}

## twelve.pep <- as.character(read.delim("/data/A_crassus/RNAseq/454_compare_diss/pep12_vs_Trinity.blt",
##                                       header = FALSE)[,2])

## contig.geno[rownames(contig.geno)%in%twelve.pep, ]


## endopeptidases <- c("comp45226_c0_seq1", "comp50970_c0_seq1", "comp51126_c0_seq1", 
##                     "comp54429_c1_seq1", "comp54639_c0_seq1", "comp55225_c0_seq1", 
##                     "comp55309_c0_seq1", "comp56603_c2_seq1", "comp57109_c0_seq1", 
##                     "comp58217_c1_seq1", "comp58550_c0_seq1", "comp58628_c0_seq1", 
##                     "comp58854_c0_seq1", "comp58859_c3_seq1", "comp58974_c1_seq1"
##                     )

## sig.endo <- endopeptidases[endopeptidases%in%sigGenes(MF.seperate[[1]])]


## contig.geno[endopeptidases, ]

## contig.geno[sig.endo, ]

## table(sig.endo%in%twelve.pep)

## both <- twelve.pep[(twelve.pep%in%sig.endo)]

## contig.geno[both, ]

BLnames <- c("454_vs_Bm.blt", 
             "454_vs_Ce.blt",
             "a_Trinity_vs_Bm.blt",
             "a_Trinity_vs_Ce.blt",
             "Trinity_vs_Bm.blt",
             "Trinity_vs_Ce.blt")

BL <- vector("list", length(BLnames))
names(BL) <- BLnames


CR <- list()
for (w in names(BL)){
    blpath <- paste("/data/A_crassus/RNAseq/454_compare_diss/", w, sep = "")
    if (grepl("Bm", w)){dbpath = "/data/db/blastdb/uniref100/uniref100_bm.fasta"}
    if (grepl("Ce", w)){dbpath = "/data/db/blastdb/wormpep220/wormpep220.fasta"}
    command <- paste("/tools/blast_mask_myStrange_Fasta.pl -d",
                     dbpath, blpath ," -co 'L' 2>&1 1>/dev/null")
    CR[[w]] <- read.delim(pipe(command),
                          header=FALSE)
}

CR <- melt(CR, id.vars=c("V1", "V2", "V3","V4"))

tapply(CR$V3, CR$L1, sum) /( tapply(CR$V4, CR$L1, sum)  + tapply(CR$V3, CR$L1, sum) )



allFFF <- names(read.ESeq("/data/db/blastdb/A_crassus/454_assembly.fasta"))
allFFF <- gsub(" .*$", "", allFFF)

FFFvsTri <- read.delim("/data/A_crassus/RNAseq/454_compare_diss/454_vs_Trinity/454_vs_Trinity.blt", header = FALSE)
names(FFFvsTri)[1:2] <- c("FFF", "Tri")

allTri <- names(read.ESeq("/data/db/blastdb/A_crassus/Trinity.fasta"))
allTri <- gsub(" .*$", "", allTri)

TrivsFFF <- read.delim("/data/A_crassus/RNAseq/454_compare_diss/454_vs_Trinity/Trinity_vs_454.blt", header = FALSE)
names(TrivsFFF)[1:2] <- c("Tri", "FFF")

presence <- merge(FFFvsTri[, 1:2], TrivsFFF [, 1:2], all = TRUE )
##presence$Tri <- as.character(presence$Tri)

presence <- merge(presence, as.data.frame(allTri), by.x = "Tri",
                  by.y = "allTri", all = TRUE)

presence <- merge(presence, as.data.frame(allFFF), by.x = "FFF",
                  by.y = "allFFF", all = TRUE)

presence.num <- as.data.frame(apply(presence, 2,
                                    function (x) (as.numeric(!is.na(x)))))
vennDiagram(presence.num)

BLT.FFFvsTri <- read.delim("/data/A_crassus/RNAseq/454_compare_diss/454_vs_Trinity/454_Trinity.blmask", header = FALSE)

BLT.TrivsFFF <- read.delim("/data/A_crassus/RNAseq/454_compare_diss/454_vs_Trinity/Trinity_vs_454.blmask", header = FALSE)



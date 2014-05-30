library(rtracklayer)
library(biomaRt)

if(!exists("T.e")){
  source("src/3_edgeR.R")
}

gff <- import.gff3("/data/A_crassus/RNAseq/protein_prediction/best_candidates.eclipsed_orfs_removed.gff3",
                   asRangedData=FALSE)

contig2gene <- subset(as.data.frame(gff), type%in%"mRNA")
contig2gene <- contig2gene[, c("seqnames", "ID")]
names(contig2gene) <- c("transcript", "protein")
contig2gene$gene <- gsub("_seq\\d+", "", contig2gene$transcript)

#### NR annotation
nuc_annot <- read.delim("/data/A_crassus/RNAseq/BLAST/nuc_vs_nr/Trinity.fasta.vs.nr.bltTaxAnnot",
                        comment.char = "#", header=FALSE)
names(nuc_annot) <- c("transcript", "hit", "evalue", "taxid", "nuc_nr_annot")
nuc_annot <- nuc_annot[!duplicated(nuc_annot$transcript),]

## filtering for taxon and low coverage
nuc_annot <- nuc_annot[nuc_annot$transcript%in%good.Tax.transripts, ]
nuc_annot <- nuc_annot[nuc_annot$transcript%in%rownames(T.e), ]
##############################

#### SwissProt annotation
sprot_annot <- read.delim("/data/A_crassus/RNAseq/BLAST/pep_vs_sprot/best_candidates.blt_best_annot",
                          comment.char = "#", header=FALSE)
names(sprot_annot) <- c("protein", "hit", "evalue",
                        "qstart", "qend", "sstart", "send",
                        "prot_sprot_annot")
sprot_annot <- sprot_annot[!duplicated(sprot_annot$protein), ]
sprot_annot$id <- gsub(".*\\|(.*)\\|.*", "\\1", sprot_annot$hit)

sprot_annot <- merge(sprot_annot, contig2gene)

## filtering for taxon and low coverage
sprot_annot <- sprot_annot[sprot_annot$transcript%in%good.Tax.transripts, ]
sprot_annot <- sprot_annot[sprot_annot$transcript%in%rownames(T.e), ]
##############################

### merge to one blast annotation #### 
b_annot <- merge(nuc_annot[, c("transcript", "nuc_nr_annot")],
                 sprot_annot[, c("transcript", "prot_sprot_annot", "protein", "gene")],
                 by= "transcript", all = TRUE)

#### IPR/GO annotation

## GO via sprot
if(!file.exists("/data/A_crassus/RNAseq/annotation/Uniprot_GO_EC.tsv")){
    uniProt <- useMart("unimart",dataset="uniprot",host="www.ebi.ac.uk",
                       path="/uniprot/biomart/martservice")
    all.uni.names <- getBM(attributes =c("accession", "name",
                               "protein_name", "go_id"),
                           filter="accession",
                           values=sprot_annot$id, mart=uniProt)
    write.csv(all.uni.names,
              "/data/A_crassus/RNAseq/annotation/Uniprot_GO_EC.tsv",
              row.names = FALSE)
} else {
    all.uni.names <-
        read.csv("/data/A_crassus/RNAseq/annotation/Uniprot_GO_EC.tsv")
}

sprot_GO_annot<- merge(sprot_annot, all.uni.names,
                       by.x = "id", by.y = "accession",
                       all.x = TRUE)
             
### ## GO via IPR
split.IPR.transcripts <- function (ipr){
  ## split duplicate proteins "TCONS_00010379|TCONS_00010380" over two rows
  ipr$transcript <- gsub("\\|.*", "", ipr$V1)
  ipr.alt <- ipr[grepl("\\|", ipr$V1), ]
  ipr.alt$transcript <- gsub(".*\\|", "", ipr.alt$V1)
  i <- rbind(ipr, ipr.alt)
  i$V1 <- NULL
  return(i)
}

read.ipr.frame <- function (path){
    T.lines <- readLines(path)
    T <- lapply(T.lines, function(x) unlist(strsplit(x, "\\t")))
    eleven <- as.data.frame(do.call(rbind, T[unlist(lapply(T, length))==11]))
    eleven <- cbind(eleven, NA, NA, NA, NA)
    fourteen <- as.data.frame(do.call(rbind, T[unlist(lapply(T, length))==14]))
    fourteen <- cbind(fourteen, NA)
    fiveteen <- as.data.frame(do.call(rbind, T[unlist(lapply(T, length))==15]))
    names(eleven) <- names(fiveteen)
    names(fourteen) <- names(fiveteen)
    IPR <- rbind(eleven, fourteen, fiveteen)
    split.IPR.transcripts(IPR)
}

IPR <- read.ipr.frame("/data/A_crassus/RNAseq/IPR/best_candidates.eclipsed_orfs_removed.IPR.tsv")

IPR$transcript <- gsub(".* (comp\\d+_c\\d+_seq\\d+)\\:.*", "\\1",
                       IPR$transcript)
IPR <- merge(IPR, contig2gene)

## ## filtering for taxon and low coverage
IPR <- IPR[IPR$transcript%in%good.Tax.transripts, ]
IPR <- IPR[IPR$transcript%in%rownames(T.e), ]

get.GO.df <- function (ipr, sum.col, name){
    GO.list <- by(ipr, sum.col, function (x){
        unlist(strsplit(unique(as.character(x$V14)), "\\|"))})
    GO.list <- GO.list[lapply(GO.list, length)>0]
    rep.go <- lapply(GO.list, length)
    GO.df <- as.data.frame(unlist(GO.list))
    GO.df$gene <- rep(names(GO.list), times = rep.go)
    names(GO.df) <- c("go_id", name)
    GO.df <- GO.df[!is.na(GO.df$go_id), ]
    GO.df$go_id<- as.character(GO.df$go_id)    
    rownames(GO.df) <- NULL
    return(GO.df)
}

GO.df.gene <- get.GO.df(IPR, IPR$gene, "gene")
GO.df.gene$method <- "IPR"

## GO.df.gene <- rbind(GO.df.gene,
##                     cbind(sprot_GO_annot[,c("gene", "go_id")],
##                           method = "blast"))

GO.df.transcript <- get.GO.df(IPR, IPR$transcript, "transcript")
GO.df.transcript$method <- "IPR"

## GO.df.transcript <- rbind(GO.df.transcript,
##                      cbind(sprot_GO_annot[,c("transcript", "go_id")],
##                            method = "blast"))

gene.2.GO <- by(GO.df.gene, GO.df.gene$gene,
                function (x) c(unique(as.character(x$go_id))))

## GO.2.gene <- by(GO.df.gene, GO.df.gene$go,
##                 function (x) c(as.character(x$gene)))

transcript.2.GO <- by(GO.df.transcript, GO.df.transcript$transcript,
                      function (x) c(unique(as.character(x$go_id))))

## GO.2.transcript <- by(GO.df.transcript, GO.df.transcript$go,
##                 function (x) c(as.character(x$gene)))


## SignalP via iprscan
SS.genes <- unique(IPR[IPR$V4%in%"SignalP_EUK" & IPR$V5%in%"SignalP-noTM", "gene"])
TM.genes <- unique(IPR[IPR$V4%in%"SignalP_EUK" & IPR$V5%in%"SignalP-TM", "gene"])

SS.trans <- unique(IPR[IPR$V4%in%"SignalP_EUK" & IPR$V5%in%"SignalP-noTM", "transcript"])
TM.trans <- unique(IPR[IPR$V4%in%"SignalP_EUK" & IPR$V5%in%"SignalP-TM", "transcript"])

TMHMM.genes <- unique(IPR[IPR$V4%in%"TMHMM", "gene"])
TMHMM.trans <- unique(IPR[IPR$V4%in%"TMHMM", "transcript"])


library(rtracklayer)

if(!exists("good.Tax.transripts")){
  source("/home/ele/thesis/experimental_infection/basic/taxon_screen.R")
}


gff <- import.gff3("/data/RNAseq/protein_prediction/best_candidates.eclipsed_orfs_removed.gff3",
                   asRangedData=FALSE)

contig2gene <- subset(as.data.frame(gff), type%in%"mRNA")
contig2gene <- contig2gene[, c("seqnames", "ID")]
names(contig2gene) <- c("transcript", "protein")
contig2gene$gene <- gsub("_seq\\d+", "", contig2gene$transcript)

nuc_annot <- read.delim("/data/RNAseq/BLAST_nuc_vs_nr/Trinity.fasta.vs.nr.bltTaxAnnot",
                        comment.char = "#", header=FALSE)
names(nuc_annot) <- c("gene", "hit", "evalue", "taxid", "nuc_nr_annot")
nuc_annot$name <- gsub("_seq\\d+", "", nuc_annot$gene)
nuc_annot <- nuc_annot[!duplicated(nuc_annot$gene),]

sprot_annot <- read.delim("/data/RNAseq/BLAST_pep_vs_sprot/best_candidates.blt_best_annot",
                          comment.char = "#", header=FALSE)
names(sprot_annot) <- c("protein", "hit", "evalue",
                        "qstart", "qend", "sstart", "send",
                        "prot_sprot_annot")
sprot_annot <- sprot_annot[!duplicated(sprot_annot$protein), ]


b_annot <- merge(contig2gene, nuc_annot[,c("gene", "nuc_nr_annot")],
                 by.x = "transcript", by.y = "gene", all.x = TRUE)

b_annot <- merge(b_annot, sprot_annot[,c("protein", "prot_sprot_annot")],
                 by.x = "protein", by.y = "protein", all.x = TRUE)

b_annot$transcript <- as.character(b_annot$transcript)
head(gsub("[\\d+|\\w+|\\|]", "", b_annot$prot_sprot_annot))

b_annot <- b_annot[b_annot$transcript%in%good.Tax.transripts,]

IPR <- read.delim("/data/RNAseq/IPR/best_candidates.eclipsed_orfs_removed.IPR.tsv",
                  header=FALSE)

IPR$transcript <- gsub(".* (comp\\d+_c\\d+_seq\\d+)\\:.*", "\\1", IPR$V1)
IPR$gene <- gsub(".* (comp\\d+_c\\d+)_seq\\d+\\:.*", "\\1", IPR$V1)
IPR$protein <- gsub("(^m\\.\\d+).*", "\\1", IPR$V1)

get.GO.df <- function (ipr, sum.col){
  GO.list <- by(ipr, sum.col, function (x){
    unlist(strsplit(unique(as.character(x$V14)), "\\|"))})
  GO.list <- GO.list[lapply(GO.list, length)>0]
  rep.go <- lapply(GO.list, length)
  GO.df <- as.data.frame(unlist(GO.list))
  GO.df$gene <- rep(names(GO.list), times = rep.go)
  names(GO.df) <- c("go", "gene")
  GO.df$go <- as.character(GO.df$go)
  rownames(GO.df) <- NULL
  return(GO.df)
}

GO.df.gene <- get.GO.df(IPR, IPR$gene)
GO.df.transcript <- get.GO.df(IPR, IPR$transcript)

gene.2.GO <- by(GO.df.gene, GO.df.gene$gene,
                function (x) c(as.character(x$go)))

GO.2.gene <- by(GO.df.gene, GO.df.gene$go,
                function (x) c(as.character(x$gene)))

transcript.2.GO <- by(GO.df.transcript, GO.df.transcript$gene,
                      function (x) c(as.character(x$go)))

GO.2.transcript <- by(GO.df.transcript, GO.df.transcript$go,
                function (x) c(as.character(x$gene)))

library(rtracklayer)

if(!exists("good.Tax.transripts")){
  source("basic/taxon_screen.R")
}

gff <- import.gff3("/data/A_crassus/RNAseq/protein_prediction/best_candidates.eclipsed_orfs_removed.gff3",
                   asRangedData=FALSE)

contig2gene <- subset(as.data.frame(gff), type%in%"mRNA")
contig2gene <- contig2gene[, c("seqnames", "ID")]
names(contig2gene) <- c("transcript", "protein")
contig2gene$gene <- gsub("_seq\\d+", "", contig2gene$transcript)

nuc_annot <- read.delim("/data/A_crassus/RNAseq/BLAST/nuc_vs_nr/Trinity.fasta.vs.nr.bltTaxAnnot",
                        comment.char = "#", header=FALSE)
names(nuc_annot) <- c("gene", "hit", "evalue", "taxid", "nuc_nr_annot")
nuc_annot$name <- gsub("_seq\\d+", "", nuc_annot$gene)
nuc_annot <- nuc_annot[!duplicated(nuc_annot$gene),]

sprot_annot <- read.delim("/data/A_crassus/RNAseq/BLAST/pep_vs_sprot/best_candidates.blt_best_annot",
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


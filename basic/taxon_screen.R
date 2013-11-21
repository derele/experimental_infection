## Read blast reports for annotation and Taxon screening 
nuc_Tax <- read.csv("/data/A_crassus/RNAseq/BLAST/nuc_vs_nr/Trinity.fasta.vs.nr.ALLTAX",
                    sep = ",", as.is=TRUE, header=FALSE)

names(nuc_Tax) <- c("transcript", "taxid",
                    "family", "phylum",
                    "kingdom" ,"superkingdom")
nuc_Tax <- nuc_Tax[!duplicated(nuc_Tax$transcript),]
nuc_Tax$gene <- gsub("_seq\\d+", "", nuc_Tax$transcript)


nuc_Tax_N <- read.csv("/data/A_crassus/RNAseq/BLAST/nuc_vs_nt/Trinity.fasta_vs_nt.ALLTAX",                      sep = ",", as.is=TRUE, header=FALSE)

names(nuc_Tax_N) <- c("transcript", "taxid",
                      "family", "phylum",
                      "kingdom" ,"superkingdom")

nuc_Tax_N<- nuc_Tax_N[!duplicated(nuc_Tax_N$transcript),]
nuc_Tax_N$gene <- gsub("_seq\\d+", "", nuc_Tax_N$transcript)


## Think about how to weight nt vs. nr (bitscore ?)

good.Tax.genes <- unique(nuc_Tax[!nuc_Tax$phylum%in%"Chordata", "gene"])
good.Tax.transripts <- unique(nuc_Tax[!nuc_Tax$phylum%in%"Chordata",
                                      "transcript"])



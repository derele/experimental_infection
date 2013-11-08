



## Read blast reports for annotation and Taxon screening 
nuc_Tax <- read.csv("/data/RNAseq/BLAST_nuc_vs_nr/Trinity.fasta.vs.nr.ALLTAX",
                    sep = ",", as.is=TRUE)

names(nuc_Tax) <- c("transcript", "taxid",
                    "family", "phylum",
                    "kingdom" ,"superkingdom")
nuc_Tax$gene <- gsub("_seq\\d+", "", nuc_Tax$transcript)

good.Tax.genes <- unique(nuc_Tax[!nuc_Tax$phylum%in%"Chordata", "gene"])
good.Tax.transripts <- unique(nuc_Tax[!nuc_Tax$phylum%in%"Chordata",
                                      "transcript"])

## prot_Tax

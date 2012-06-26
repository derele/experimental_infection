
write.by.annotation <- function (annotation, path){
  selected.contigs <- contig.df[grepl(annotation, contig.df$Bm.annot) &
                            contig.df$category%in%"MN", c("contig", "mstrand", "seq")]
  selected.fasta <- apply(selected.contigs, 1,
                      function (x) if (x["mstrand"]){revcom(x["seq"])}
                      else {as.character(x["seq"])})
  names(selected.fasta) <- selected.contigs[["contig"]]
  write.sequence(selected.fasta, path)
}

write.by.contig <- function (contig, path){
  selected.contigs <- contig.df[contig.df$contig%in%contig &
                            contig.df$category%in%"MN", c("contig", "mstrand", "seq")]
  selected.fasta <- apply(selected.contigs, 1,
                      function (x) if (x["mstrand"]){revcom(x["seq"])}
                      else {as.character(x["seq"])})
  names(selected.fasta) <- selected.contigs[["contig"]]
  write.sequence(selected.fasta, path)
}





write.by.contig("Contig1879",
                "/home/ele/thesis/experimental_infection/q_pcr/Cox1.fasta")

write.by.contig("Contig6759",
                "/home/ele/thesis/experimental_infection/q_pcr/Cox2.fasta")

write.by.annotation("Cytochrome.*subunit 3",
                    "/home/ele/thesis/experimental_infection/q_pcr/Cox3.fasta")

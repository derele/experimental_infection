fastq <- list.files("/data/A_crassus/RNAseq/renamed/",
                    pattern = "left", full.names = TRUE)


raw.reads <- list()
for(file in fastq){
    raw.reads[file] <- readLines(pipe(paste("grep @", file, "| wc -l")))
}


raw.reads <- unlist(raw.reads)
names(raw.reads) <- gsub("/data/A_crassus/RNAseq/renamed//(\\w{2}_.*)_left.fastq",
                         "\\1", names(raw.reads))

if (!exists("sex.conds")){
    source("exp_test/edgeR.R")
}

## Reading the Rsem couts for transcripts
ntrans <- read.delim("/data/A_crassus/RNAseq/rsem_trinity/Trinity_trans.counts.matrix",
                     as.is = TRUE)
ntrans.matrix <- as.matrix(ntrans[,-1])
rownames(ntrans.matrix) <- ntrans[,1]
colnames(ntrans.matrix) <-  gsub("\\.isoforms\\.results", "", colnames(ntrans.matrix))


## mapping to all genes
mapped.reads <- round(colSums(ntrans.matrix))
analysed.reads <- round(colSums(ntrans.matrix[rownames(T.e), ]))

host.species <- ifelse(eel.conds%in%"AA", "An. anguilla", "An. japonica")

n.womrs.preped <- c(AA_R11M = 14, AA_R16M = 4, AA_R18F = 1, AA_R28F = 1, AA_R2M = 4, AA_R8F = 1, 
                    AA_T12F = 1, AA_T20F = 1, AA_T24M = 3, AA_T3M = 4, AA_T42M = 1, AA_T45F = 1, 
                    AJ_R1F = 1, AJ_R1M = 1, AJ_R3F = 1, AJ_R3M = 2, AJ_R5F = 1, AJ_R5M = 1,
                    AJ_T19M = 7, AJ_T20M = 8, AJ_T25M = 5, AJ_T26F = 1, AJ_T5F = 1, AJ_T8F = 1)

worm.population <- c(AA_R11M = "Europe (R)", AA_R16M = "Europe (R)", AA_R18F = "Europe (R)",
                     AA_R28F = "Europe (R)", AA_R2M = "Europe (B)", AA_R8F = "Europe (B)", 
                     AA_T12F = "Taiwan (K)", AA_T20F = "Taiwan (K)", AA_T24M = "Taiwan (K)",
                     AA_T3M = "Taiwan (Y)", AA_T42M = "Taiwan (Y)", AA_T45F = "Taiwan (Y)", 
                     AJ_R1F = "Europe (R)", AJ_R1M = "Europe (R)", AJ_R3F = "Europe (R)",
                     AJ_R3M = "Europe (R)", AJ_R5F = "Europe (B)", AJ_R5M = "Europe (B)",
                     AJ_T19M = "Taiwan (Y)", AJ_T20M = "Taiwan (Y)", AJ_T25M = "Taiwan (Y)",
                     AJ_T26F = "Taiwan (Y)", AJ_T5F = "Taiwan (K)", AJ_T8F = "Taiwan (Y)")


table1 <- cbind(exp.host.species= host.species, worm.sex=as.character(sex.conds), worm.population,
                n.womrs.preped, raw.reads, mapped.reads, analysed.reads)

library(xtable)

print(xtable(table1), type="html", file="/home/ele/Dropbox/Ac_trans_div/tables/Table1.html")


if(!exists("dapc1")){
    source("exp_test/edgeR.R")
}

diff.SNPs <- rownames(GT[dapc1$pca.loadings[,1]>0.006, ])
#diff.SNPs <- rownames(GT[dapc1$pca.loadings[,1]>0.01, ])

diff.SNP.genes <- gsub("_seq\\d+\\.\\d+", "", diff.SNPs)

diff.SNP.genes <- diff.SNP.genes[diff.SNP.genes%in%rownames(G.e)]


## No association
SNP.geno.pop <- table(rownames(G.e)%in%diff.SNP.genes,
                      rownames(G.e)%in%g.contigs.l[[3]])
fisher.test(SNP.geno.pop)


SNP.geno.eel <- table(rownames(G.e)%in%diff.SNP.genes,
                      rownames(G.e)%in%g.contigs.l[[2]])
fisher.test(SNP.geno.eel)

SNP.geno.sex <- table(rownames(G.e)%in%diff.SNP.genes,
                      rownames(G.e)%in%g.contigs.l[[1]])
fisher.test(SNP.geno.sex)


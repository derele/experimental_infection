

perl -ne '/DP=(\d*)/; print($1,"\n")'  genotypes.vcf
awk '{if ($6 > 20) {print $0}}' genotypes.vcf | head -n 20 | column -t | most

awk '{if ($6 > 20) {print $0}}' genotypes.vcf > genotypes_above30.vcf
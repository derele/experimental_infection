#!/bin/bash

cd /data/RNAseq/;
samtools mpileup -uf /data/RNAseq/rsem_trinity/*.sorted.bam | bcftools view -gcv -

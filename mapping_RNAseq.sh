#!/bin/bash

# ------------- system commands used by this script --------------------
BWA=/home/ele/tools/bwa-0.5.9/bwa;
SAMT=/home/ele/tools/samtools/samtools

# ------------- the files that should be mapped     --------------------

find .. -name "*.solfastq.gz" | parallel "gunzip -c {} | /home/ele/tools/fq_all2std.pl sol2std > {}.fastq";

find .. -name "*.fastq" | parallel "$BWA aln db_assembly/fullest_assembly.fasta {} > {}.sai";

function exitOnError
  {
  local EXITCODE="$?"
  if [ $EXITCODE -ne "0" ]; then
    echo -e "\nERROR: '$1' failed with error code: '$EXITCODE'\n"
    exit 1
  fi
  }

SA=($(find .. -name "*.fastq.sai" | tr ' ' '\n' | sort -r));

FA=($(find .. -name "*.fastq" | tr ' ' '\n' | sort -r));

for ((i=0; i< ${#SA[*]}; i=i+4));
do 

echo "$BWA sampe db_assembly/fullest_assembly.fasta ${SA[$i]} ${SA[$i+1]} ${FA[$i]} ${FA[$i+1]} | $SAMT view -uS -t db_assembly/fullest_assembly.fasta.fai  - | $SAMT sort -  $(dirname ${SA[$i]} | awk '{sub("../", "", $1); print $1}')_1;"

$BWA sampe db_assembly/fullest_assembly.fasta ${SA[$i]} ${SA[$i+1]} ${FA[$i]} ${FA[$i+1]} | $SAMT view -uS -t db_assembly/fullest_assembly.fasta.fai  - | $SAMT sort -  $(dirname ${SA[$i]} | awk '{sub("../", "", $1); print $1}')_1 exitOnError; "pipe1"

echo "$BWA sampe db_assembly/fullest_assembly.fasta ${SA[$i+2]} ${SA[$i+3]} ${FA[$i+2]} ${FA[$i+3]} | $SAMT view -uS -t db_assembly/fullest_assembly.fasta.fai  - | $SAMT sort -  $(dirname ${SA[$i+2]} | awk '{sub("../", "", $1); print $1}')_2;"

$BWA sampe db_assembly/fullest_assembly.fasta ${SA[$i+2]} ${SA[$i+3]} ${FA[$i+2]} ${FA[$i+3]} | $SAMT view -uS -t db_assembly/fullest_assembly.fasta.fai  - | $SAMT sort -  $(dirname ${SA[$i+2]} | awk '{sub("../", "", $1); print $1}')_2; exitOnError "pipe2"

done;


## indes all the created files
find -name "*.bam" | parallel /home/ele/tools/samtools/samtools index; exitOnError "index" 

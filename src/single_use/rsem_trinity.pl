#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);
use FindBin;
use Cwd;

######################################################
## Set to base directory of the Trinity installation:
my $BASEDIR = "/home/ele/tools/trinityrnaseq_r2013-02-16/";
######################################################

#############################################################
## Run RSEM, compute abundance estimates for each replicate.
#############################################################

my @trans_rsem_files;
my @genes_rsem_files;

my @replicates = qw(AA_R11M AA_R16M AA_R18F AA_R28F AA_R2M AA_R8F AA_T12F AA_T20F AA_T24M AA_T3M AA_T42M AA_T45F AJ_R1F AJ_R1M AJ_R3F AJ_R3M AJ_R5F AJ_R5M AJ_T19M AJ_T20M AJ_T25M AJ_T26F AJ_T5F AJ_T8F);

foreach my $replicate (@replicates) {
  my $left_fq_file = "../renamed/".$replicate."_left.fastq";
  my $right_fq_file = "../renamed/".$replicate."_right.fastq";
        
  ## run RSEM
  ##         my $seqType = $PARAMS{"--seqType"} or die "Error, --seqType not specified";
  my $cmd = "$BASEDIR/util/RSEM_util/run_RSEM_align_n_estimate.pl "
    . " --transcripts /data/RNAseq/Trinity/Trinity.fasta "
      . " --seqType fq "
        . " --prefix $replicate "
          ." --left $left_fq_file --right $right_fq_file "
            ."--thread_count 8";
  
  &process_cmd($cmd,
               "Running RSEM to perform abundance estimation on sample $replicate"
              ) unless (-s "$replicate.isoforms.results");
        
  unless (-s "Trinity.trans_lengths.txt") {
            
    ## Create the effective lengths file:
    ## just need to process one of the reps for this.  These data should be identical among all such files.

    &process_cmd("cat $replicate.isoforms.results | cut -f1,3,4 > Trinity.trans_lengths.txt",
                 "Generating Trinity transcript lengths file");
    &process_cmd("cat $replicate.genes.results | cut -f1,3,4 > Trinity.components_lengths.txt",
                 "Generating Trinity component lengths file");
  }
        
  push (@trans_rsem_files, "$replicate.isoforms.results");
  push (@genes_rsem_files, "$replicate.genes.results");
}

######################################
# Make count matrices for DE analysis
######################################

# # make iso matrix
my $cmd = "$BASEDIR/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl " . join (" " , @trans_rsem_files) . " > Trinity_trans.counts.matrix";
&process_cmd($cmd,
             "Generating transcript (isoform) count matrix file."
            ) unless (-s "Trinity_trans.counts.matrix");

# # # make the 'gene' matrix
$cmd = "$BASEDIR/util/RSEM_util/merge_RSEM_frag_counts_single_table.pl " . join (" " , @genes_rsem_files) . " > Trinity_components.counts.matrix"; 
&process_cmd($cmd,
             "Generating Trinity component (gene) count matrix file"
            ) unless (-s "Trinity_components.counts.matrix");

################################
## Run the TMM normalization:
################################

$cmd = "$BASEDIR/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix Trinity_trans.counts.matrix --lengths Trinity.trans_lengths.txt ";
&process_cmd($cmd, "Running TMM scaling normalization and generating FPKM matrix for transcripts (isoforms)");

$cmd = "$BASEDIR/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix Trinity_components.counts.matrix --lengths Trinity.components_lengths.txt ";
&process_cmd($cmd, "Running TMM scaling normalization and generating FPKM matrix for components (genes)");


# ################
# ## run EdgeR
# ################

# foreach my $target_type ("trans", "components") {
    
#     my $edgeR_dir = "edgeR_${target_type}";

#     my $cmd = "$BASEDIR/Analysis/DifferentialExpression/run_DE_analysis.pl "
#         . " --matrix Trinity_${target_type}.counts.matrix "
#         . " --method edgeR "
#         . " --samples_file $read_samples_descr_file "
#         . " --output $edgeR_dir ";
    
#     &process_cmd($cmd, "Running edgeR for $target_type") unless (-d $edgeR_dir);

    
#     chdir $edgeR_dir or die "Error, cannot cd to $edgeR_dir";
    
#     ## extract the diff. expressed transcripts.
#     $cmd = "$BASEDIR/Analysis/DifferentialExpression/analyze_diff_expr.pl "
#         . " --matrix ../Trinity_${target_type}.counts.matrix.TMM_normalized.FPKM ";
    

#     if (exists $PARAMS{"-P"}) {
#         $cmd .= " -P " . $PARAMS{"-P"};
#     }
#     if (exists $PARAMS{"-C"}) {
#         $cmd .= " -C " . $PARAMS{"-C"};
#     }
    
#     &process_cmd($cmd, "Running analysis of DE $target_type, hierarchically clustering transcripts and samples");
    
#     chdir $workdir or die "Error, cannot cd back to $workdir";
# }


# print "Done.\n";
    



exit(0);


####
sub process_cmd {
    my ($cmd, $msg) = @_;


    if ($msg) {
        print "\n\n";
        print "#################################################################\n";
        print "$msg\n";
        print "#################################################################\n";
    }
    
    print "CMD: $cmd\n";

    my $time_start = time();
    
    my $ret = system($cmd);
    my $time_end = time();

    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }

    my $number_minutes = sprintf("%.1f", ($time_end - $time_start) / 60);
    
    print "TIME: $number_minutes min. for $cmd\n";
    

    return;
}


# ####
# sub parse_sample_descriptions {
#     my ($read_samples_descr_file, $PARAMS_href) = @_;

#     my %samples_descr;
    
    
#     open (my $fh, $read_samples_descr_file) or die $!;
#     while (<$fh>) {
#         if (/^\#/) { next; }
#         unless (/\w/) { next; }
#         s/^\s+|\s+$//g;
#         chomp;
#         my @x = split(/\t/);
#         if (/=/ && scalar(@x) == 1) {
#             my ($key, $val) = split(/=/, $x[0]);
#             $PARAMS_href->{$key} = $val;
#         }
#         else {
        
            
#             my ($condition, $replicate, $reads_left, $reads_right) = @x;
            
#             ## remove gzip extension, will convert to gunzipped version later 
#             $reads_left =~ s/\.gz$//;
#             $reads_right =~ s/\.gz$// if $reads_right;
        
#             $samples_descr{$condition}->{$replicate} = [$reads_left, $reads_right];
#         }
#     }
    
#     close $fh;
    
#     return(%samples_descr);
    
    
# }

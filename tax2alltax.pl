#!/usr/bin/perl

use warnings;      # avoid d'oh! bugs
use strict;        # avoid d'oh! bugs
$|=1;
use Data::Dumper;  # easy printing, sometimes only used during coding

use Bio::Perl;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;
use Bio::LITE::Taxonomy::NCBI::Gi2taxid; # qw/new_dict/;

# new_dict (in => "/db/blastdb/taxonomy/gi_taxid_prot.dmp",
#            out => "/db/blastdb/taxonomy/gi_taxid_prot.bin");  

my $taxDB = Bio::LITE::Taxonomy::NCBI-> new (
                                             db=>"NCBI",
                                             names=>
                                             "/db/blastdb/taxonomy/names.dmp",
                                             nodes=>
                                             "/db/blastdb/taxonomy/nodes.dmp",
                                             dict=>"/db/blastdb/taxonomy/gi_taxid_prot.bin"
                                            );

while (<>) {
  next if /#/;
  chomp($_);
  my @blt = split (" ", $_);
  my $contig = $blt[0];
  my $taxid = $blt[3];
  my $family = $taxDB->get_term_at_level($taxid,"family");
  my $phylum = $taxDB->get_term_at_level($taxid,"phylum");
  my $kingdom = $taxDB->get_term_at_level($taxid,"kingdom");
  my $superkingdom = $taxDB->get_term_at_level($taxid,"superkingdom");

  print "$contig,$taxid,$family,$phylum,$kingdom,$superkingdom\n";
}


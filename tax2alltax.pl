#!/usr/bin/perl

use warnings;      # avoid d'oh! bugs
use strict;        # avoid d'oh! bugs
$|=1;
use Data::Dumper;  # easy printing, sometimes only used during coding

use Bio::Perl;
use Bio::LITE::Taxonomy;
use Bio::LITE::Taxonomy::NCBI;

my $taxDB = Bio::LITE::Taxonomy::NCBI-> new (
                                             db=>"NCBI",
                                             names=>
                                             "/data/db/taxonomy/names.dmp",
                                             nodes=>
                                             "/data/db/taxonomy/nodes.dmp",
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


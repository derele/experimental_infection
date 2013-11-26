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

my $cmd = "blast_formatter -archive ".$ARGV[0].
  " -outfmt \'6 qseqid sacc bitscore staxids stitle\'"; 

print STDERR "running $cmd \n";

my @blast = `$cmd`; 


print STDERR "obtaining family, phylum, kingdom and superkingdom\n\n";

print  "query,subject,bitscore,taxid,family,phylum,kingdom,superkingdom\n";

foreach (@blast) {
  next if /#/;
  chomp($_);
  my @blt = split (" ", $_);
  my $family = $taxDB->get_term_at_level($blt[3],"family");
  my $phylum = $taxDB->get_term_at_level($blt[3],"phylum");
  my $kingdom = $taxDB->get_term_at_level($blt[3],"kingdom");
  my $superkingdom = $taxDB->get_term_at_level($blt[3],"superkingdom");

  print "$blt[0],$blt[1],$blt[2],$blt[3],$family,$phylum,$kingdom,$superkingdom\n";
}


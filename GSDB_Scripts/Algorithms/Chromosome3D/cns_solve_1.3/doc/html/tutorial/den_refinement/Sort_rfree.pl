#!/usr/bin/perl
use warnings;
use strict;
##############################################################################
#### script to rank CNS-DEN refined pdb files based on r_free
####
#### Abhinav Kumar, Debanu Das
#### JCSG (JOINT CENTER FOR STRUCTURAL GENOMICS) www.jcsg.org
#### Oct 2011
##############################################################################
#### RUN THIS SCRIPT WITHOUT ANY EDITS AFTER THE CNS-DEN CLUSTER JOB ENDS
##############################################################################
my $file_dir = "Outputs";    ## the output pdb files are located here
##############################################################################
my %rfree;
foreach my $pdb (<$file_dir/*.pdb>) {
   my $r_rfree = `grep "^REMARK final  " $pdb`;
   chomp $r_rfree;
   $r_rfree =~ s/[a-z=_]//ig;
   my ( $r, $rfree ) = split " ", $r_rfree;
   $pdb =~ s/^.*\///;
   if ($rfree) {
      $rfree{"$pdb\t\t$r"} = $rfree;
   }
}
##############################################################################
print "\tPDB\t\t\t  R\tR_free\n";
foreach my $pdb ( sort { $rfree{$a} <=> $rfree{$b} } keys %rfree ) {
   print "$pdb\t$rfree{$pdb}\n";
}
##############################################################################


#!/usr/bin/perl -w

# Parse the refine.pdb files to extract results by DEN-assisted refinements

# Author: Michael Levitt

print STDERR "results_parse.pl <names of files to analyze>\n";

foreach $fn ( @ARGV ) {

  open ( IN, $fn ) || die "Cannot open input file '$fn' so STOP\n";
 
  while ( <IN> ) { 
    chop; 
    s/\r\n//g;       # eliminate PC ctrl-M due to email
    s/\015//g;       # eliminate Mac ctrl-M due to email

#REMARK final    r= 0.2294 free_r= 0.3600
    if ( /REMARK final/ && /free_r/ ) { $final = $_; next }

#REMARK input coordinates: 3crw-cns.pdb
    $pdb="none";
    if ( /REMARK input coordinates: ([^ ]+)\.pdb/ ) { $pdb = $1; next }

#REMARK DEN is active. gamma=1.0   kappa=0.1   wden=0

    if ( /REMARK DEN is active. gamma= *([\.\d]+) +kappa= *([\.\d]+) + wden= *([\.\d]+)/ ) { $gamma = $1; $wden = $3; next }
 
# REMARK positional (xyz) minimization steps= 0
    if ( /REMARK positional \(xyz\) minimization steps= *([\.\d]+)/ ) { $steps = $1; next }

# REMARK starting temperature= 3000  total md steps= 60 * 6  seed= 1
    if ( /REMARK starting temperature=.+seed= +(\d+)/ ) { $seed = $1; next  }

  }

  $pdb = substr($pdb,0,4);
  print "$pdb";
  if ( $steps == 0 ) { print ", no minimization\n" }
  else               { print ", minimization\n" }

  print "gamma=$gamma wden=$wden";
  if ( $wden == 0 ) { print " (i.e., no den)\n" } else { print "\n" }

  print "${seed}/refine.pdb:$final fn=$fn\n";

  close ( IN );
}

#!/usr/bin/perl -w

# sort the results tabulated by the perl script "results_parse.pl"

# Author: Michael Levitt

my ( $fn, $ARGV, $il, $line1, $line2, $line3, $path, $rfact, $free_r,
     $gamma, $wden, $seed, $rest, $key, $key0, $x, $com, $cond,
     @wdens,
     %Rfree_n, %Rfree_s1, %Rfree_s2, %Rfree_mn, %Rfree_mx, %Rfree_Rf
    );
print STDERR "results_sort.pl <names of a.dat files to analyze>\n";

foreach $fn ( @ARGV ) {

  open ( IN, $fn ) || die "Cannot open input file '$fn' so STOP\n";
 
  $il = 0;
  while ( <IN> ) { 
    chop; 
    s/\r\n//g;       # eliminate PC ctrl-M due to email
    s/\015//g;       # eliminate Mac ctrl-M due to email

#3crw, no minimization
#gamma=1.0 wden=0 (i.e., no den)
#1/refine.pdb:REMARK final    r= 0.2294 free_r= 0.3600 fn=12/refine.pdb

    $il++;
    if ( $il % 3 == 1 ) { $line1 = $_ }
    if ( $il % 3 == 2 ) { $line2 = $_ }
    if ( $il % 3 == 0 ) { 
      $line3 = $_;
      if ( /(.+)REMARK final +r= +([\.\d]+) +free_r= +([\.\d]+) +fn=(.+)$/ ) { $path = $1; $rfact = $2; $free_r = $3; $fn = $4, $path =~ s/:$// }
      $path = $fn;
      $line2 =~ s/wden=0 \(i.e., no den\)/wden=0.0/;
      $pdb="none";
      if ( $line1 =~ /([^,]+), +(.+)$/ ) { $pdb = $1; $cond = $2 }
      if ( $line2 =~ /gamma=([\.\d]+) +wden=([\.\d]+)/ ) { $gamma = $1; $wden = $2 }
      if ( $line3 =~ /(\d+)\/(.*)$/ ) { $seed = $1; $rest = $2 }
      $rest =~ s/refine\.pdb\:REMARK final//;
      printf ( "%s, gamma=%7.1f  den= %8.1f %3d %15s %s\n", $line1, $gamma, $wden, $seed, $path, $rest );

      $key = $gamma."_".$wden; 
      if ( !$Rfree_n{$key} ) { $Rfree_n{$key} = 1 } else { $Rfree_n{$key}++ }
      if ( !$Rfree_s1{$key} ) { $Rfree_s1{$key} = $free_r } else { $Rfree_s1{$key} += $free_r }
      if ( !$Rfree_s2{$key} ) { $Rfree_s2{$key} = $free_r**2 } else { $Rfree_s2{$key} += $free_r**2 }
      if ( !$Rfree_mn{$key} || $Rfree_mn{$key} > $free_r ) { $Rfree_mn{$key} = $free_r; $Rfree_pa{$key} = $path; $Rfree_Rf{$key} = $rfact; }
      if ( !$Rfree_mx{$key} || $Rfree_mx{$key} < $free_r ) { $Rfree_mx{$key} = $free_r }

      if ( $wden == 0 ) { $key0 = "wden_eq_0" } else { $key0 = "wden_ne_0" }
      if ( !$Rfree_n{$key0} ) { $Rfree_n{$key0} = 1 } else { $Rfree_n{$key0}++ }
      if ( !$Rfree_s1{$key0} ) { $Rfree_s1{$key0} = $free_r } else { $Rfree_s1{$key0} += $free_r }
      if ( !$Rfree_s2{$key0} ) { $Rfree_s2{$key0} = $free_r**2 } else { $Rfree_s2{$key0} += $free_r**2 }
      if ( !$Rfree_mn{$key0} || $Rfree_mn{$key0} > $free_r ) { $Rfree_mn{$key0} = $free_r; $Rfree_pa{$key0} = $path; $Rfree_Rf{$key0} = $rfact; $Rfree_key{$key0} = $key }
      if ( !$Rfree_mx{$key0} || $Rfree_mx{$key0} < $free_r ) { $Rfree_mx{$key0} = $free_r }
      
    }
  }
  close ( IN );

#         BEST: 3crw         0.0    100.0    10   0.3386   0.00463    0.3297  0.09360    0.3455   45/refine.pdb
  print "#BEST: PDB        GAMMA     WDEN  SEED  Rfr-avg    Rfr-SD   Rfr-min  Rfa-Rfr   Rfr-max   PATH/FILE\n";

  foreach $x ( sort ( keys ( %Rfree_s1 ) ) ) { 
    if ( $x =~ /wden/ ) { 
      $key=$Rfree_key{$x};
      ($gamma,$wden) = split (/_/,$key);
      $avg = -1; $std = -1;
      if ( $Rfree_n{$key} > 0 ) { $avg = $Rfree_s1{$key}/$Rfree_n{$key}; $std = sqrt ( abs ( $Rfree_s2{$key}/$Rfree_n{$key} - $avg*$avg ) ) }
      printf ( "BESTA: %-10s %5.1f %8.1f %5d %8.4f  %8.5f  %8.4f %8.5f  %8.4f   %s  %s\n", $pdb, $gamma, $wden, $Rfree_n{$key}, $avg, $std, $Rfree_mn{$key}, $Rfree_mn{$key}-$Rfree_Rf{$key}, $Rfree_mx{$key}, $Rfree_pa{$key}, $x );
      $com = "cp -p $Rfree_pa{$key} ${x}_refine.pdb";
      print "RUN COMMAND: $com\n";
      system ( $com );

# Copy the hkl files

      $filename=$Rfree_pa{$key}; $filename=~ s/\.pdb//;
      $com = "cp -p $filename.hkl ${x}_refine.hkl";
      print "RUN COMMAND: $com\n";
      system ( $com );

    }
    else { 
      $avg = -1; $std = -1;
      if ( $Rfree_n{$x} > 0 ) { $avg = $Rfree_s1{$x}/$Rfree_n{$x}; $std = sqrt ( abs ( $Rfree_s2{$x}/$Rfree_n{$x} - $avg*$avg ) ) }
      ($gamma,$wden) = split (/_/,$x);
      printf ( "BESTB: %-10s %5.1f %8.1f %5d %8.4f  %8.5f  %8.4f %8.5f  %8.4f   %s\n", $pdb, $gamma, $wden, $Rfree_n{$x}, $avg, $std, $Rfree_mn{$x}, $Rfree_mn{$x}-$Rfree_Rf{$x}, $Rfree_mx{$x}, $Rfree_pa{$x} );
    }
  }

  @wdens = (  "3.0", "10.0", "30.0", "100.0", "300.0" );
  for ( $gamma = 0; $gamma <= 1; $gamma += 0.2 ) { 
    foreach $wden ( @wdens ) { 

      if ( $wden == 0 ) { $x = "1.0_0.0" }
      else              { $x = sprintf ( "%3.1f_%s", $gamma, $wden ) }

#      print "TEST: '$gamma' '$wden' '$x'\n";
      $avg = -1; $std = -1;
#      if ( $Rfree_n{$x} > 0 ) { $avg = $Rfree_s1{$x}/$Rfree_n{$x}; $std = sqrt ( abs ( $Rfree_s2{$x}/$Rfree_n{$x} - $avg*$avg ) ) }
      if ( !$Rfree_mn{$x} ) {printf ( "DATA: %5.1f %8.1f  %8.4f\n", $gamma, $wden, -1 );} 
      else { printf ( "DATA: %5.1f %8.1f  %8.4f\n", $gamma, $wden, $Rfree_mn{$x} ); }
    }
    print "DATA:\n";
  }

}

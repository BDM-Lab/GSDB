#!/usr/bin/perl -w

# Parse the pdb files to extract energies

print STDERR "sort.pl <names of pdb files to analyze>\n";

$key=0;

foreach $fn ( @ARGV ) {

  open ( IN, $fn ) || die "Cannot open pdb file '$fn' so STOP\n";

  $final=0;
 
  while ( <IN> ) { 
    chop; 
    s/\r\n//g;       # eliminate PC ctrl-M due to email
    s/\015//g;       # eliminate Mac ctrl-M due to email

#REMARK coordinates from molecular dynamics Etotal= 37076.7 Enoe= 32621.2
    if ( /REMARK coordinates/ && /Etotal/ ) { $final = $_; next }
# only need to read the first line
    last
  }
  
  if ( !$final == 0) {
    $key=$key+1;
    
    $filename{$key}=$fn;
    
    if ( $final =~ /REMARK coordinates from molecular dynamics Etotal=[ ]+([\.\d]+) +Enoe=[ ]+([\.\d]+)/ ) { $Etotal{$key} = $1; $Enoe{$key} = $2; };

      }
  close ( IN );
}
# create directory with name "sorted"

  mkdir "sorted";
  
  $n=0;
  foreach $x ( sort { $Etotal{$a} <=> $Etotal{$b} } ( keys ( %Etotal ) ) ) { 
  
     $n=$n+1;
     print " $n $x $Etotal{$x}  $Enoe{$x} \n";
     
     $line="cp -p $filename{$x} sorted/sort_$n\.pdb";
#     print " $line \n";
     system ( $line ) ;
     
     $filename_root=$filename{$x}; $filename_root=~ s/\.pdb//;

     $line="cp -p $filename_root.dist sorted/sort_$n\.dist";
#     print " $line \n";
     system ( $line ) ;
  }

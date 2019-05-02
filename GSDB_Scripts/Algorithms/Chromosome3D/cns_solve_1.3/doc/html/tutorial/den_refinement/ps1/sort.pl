#!/usr/bin/perl -w

# Parse the pdb files to R values

print STDERR "sort.pl <names of pdb files to analyze>\n";

$key=0;

foreach $fn ( @ARGV ) {

  open ( IN, $fn ) || die "Cannot open pdb file '$fn' so STOP\n";

  $final=0;
 
  while ( <IN> ) { 
    chop; 
    s/\r\n//g;       # eliminate PC ctrl-M due to email
    s/\015//g;       # eliminate Mac ctrl-M due to email

#REMARK final    r= 0.2646 free_r= 0.2892
    if ( /REMARK final/ && /free_r/ ) { $final = $_; next }
  }
  
  if ( !$final == 0) {
    $key=$key+1;
    
    $filename{$key}=$fn;
    
    if ( $final =~ /REMARK final    r=[ ]+([\.\d]+) +free_r=[ ]+([\.\d]+)/ ) { $r{$key} = $1; $free_r{$key} = $2; };

      }
  close ( IN );
}
# create directory with name "sorted"

  mkdir "sorted";
  
  $n=0;
  foreach $x ( sort { $free_r{$a} <=> $free_r{$b} } ( keys ( %free_r ) ) ) { 
  
     $n=$n+1;
     print " Filename= $filename{$x} Rfree= $free_r{$x}  R= $r{$x} \n";
     
     $line="cp -p $filename{$x} sorted/sort_$n\.pdb";
#     print " $line \n";
     system ( $line ) ;
     
     $filename_root=$filename{$x}; $filename_root=~ s/\.pdb//;

     $line="cp -p $filename_root.hkl sorted/sort_$n\.hkl";
#     print " $line \n";
     system ( $line ) ;
  }

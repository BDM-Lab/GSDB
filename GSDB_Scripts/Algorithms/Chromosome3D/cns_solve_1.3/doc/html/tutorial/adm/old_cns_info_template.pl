#! /usr/bin/perl -w

##############################################################################
# Read in the keyword patterns from common.idb.
##############################################################################

die "Error: environment variable CNS_SOLVE not defined.\n"
  if (not exists $ENV{CNS_SOLVE});

$common_idb = "$ENV{CNS_SOLVE}/inputs/info/common.idb";

open COMIDB, "< $common_idb"
  or die "Error: Cannot open file $common_idb\n";

%common = ();
$index = 1;

while (<COMIDB>) {
  if (/^vars:\s+(.*)/) {
    foreach $p (split /\s+/, $1) {
      $common{$p} = $index++;
    }
  }
}

close COMIDB;


##############################################################################
# Read in the keywords from task file.
# Ignore keywords matched by pattern in common.idb.
##############################################################################

sub numerically { $a <=> $b; }

%list = ();
$index = 1;

while (<>) {
  if (/\Q{===>}\E\s*(.*)/) {
    foreach $def (split /;\s+/, $1) {
      $keyw = lc $def;
      $keyw =~ s/=.*//;
      $matchcomm = $index++;
      foreach $p (keys %common) {
        if ($keyw =~ /^$p$/) {
          $matchcomm = $p;
          last;
        }
      }
      $keyp = $keyw;
      $keyp =~ s/\.\d+\./\.\\d\+\./g;
      $keyp =~ s/\.\d+$/\.\\d\+/;
      $keyp =~ s/_\d+_/_\\d\+_/g;
      $keyp =~ s/_\d+$/_\\d\+/g;
      if (not exists $list{$keyp}) {
        $list{$keyp} = $matchcomm;
      }
      else {
        if ($list{$keyp} !~ /^\d+$/) {
          if ($list{$keyp} ne $matchcomm) {
            die "Error: ambiguous patterns  for \"$keyw\"\n";
          }
        }
      }
    }
  }
}


##############################################################################
# Print parameter names matched in common.idb.
##############################################################################

print "#" x 79, "\n";
print "# Parameter names matched in common.idb:\n";

%incommon = ();

foreach $keyp (sort keys %list) {
  if ($list{$keyp} !~ /^\d+$/ ) {
    $incommon{$common{$list{$keyp}}} = $keyp;
  }
}

foreach $index (sort numerically keys %incommon) {
  print "#     $incommon{$index}\n";
}

print "#" x 79, "\n";
print "#\n";


##############################################################################
# Print new keyword patterns.
##############################################################################

%new = ();

foreach $keyp (sort keys %list) {
  if ($list{$keyp} =~ /^\d+$/ ) {
    $new{$list{$keyp}} = $keyp;
  }
}

foreach $index (sort numerically keys %new) {
  print "vars: $new{$index}\n";
  print "info:\n";
  print "#\n";
}

__END__

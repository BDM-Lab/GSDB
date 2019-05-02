#! /usr/bin/perl -w

die "Error: environment variable CNS_SOLVE not defined.\n"
  if (not exists $ENV{CNS_SOLVE});

$common_idb = "$ENV{CNS_SOLVE}/inputs/info/common.idb";


sub array_lookup {
  my $array_ref = shift;
  my $value     = shift;
  foreach $v (@{$array_ref}) {
    return 1 if ($value eq $v);
  }
  return 0;
}


##############################################################################
# Read in the keyword patterns from common.idb.
##############################################################################

open COMIDB, "< $common_idb"
  or die "Error: Cannot open file $common_idb\n";

@common = ();

while (<COMIDB>) {
  if (/^vars:\s+(.*)/) {
    foreach $p (split /\s+/, $1) {
      push @common, [$p, []];
    }
  }
}

close COMIDB;


##############################################################################
# Read in the keywords from task file.
# Ignore keywords matched by pattern in common.idb.
##############################################################################

@list = ();
@comments = ();
$mode = "";

while (<>) {
  chomp;
  if (! $mode) {
    if (/^\{=+\s\S.*\s=+\}$/) {
      @comments = ();
    }
    elsif (/^\{([*+])\s(.*)/) {
      $mode = $1;
      $text = $2;
      if ($text =~ /(.*)\s\Q$mode\E\}\s*$/) {
        $text = $1;
        $mode = "";
      }
      $text =~ s/\s+$//;
      push @comments, $text;
    }
    elsif (/\Q{===>}\E\s*(.*)/) {
      foreach $def (split /;\s+/, $1) {
        $keyw = lc $def;
        $keyw =~ s/=.*//;
        $keyp = $keyw;
        $keyp =~ s/\.\d+\./\.\\d\+\./g;
        $keyp =~ s/\.\d+$/\.\\d\+/;
        $keyp =~ s/_\d+_/_\\d\+_/g;
        $keyp =~ s/_\d+$/_\\d\+/g;
        $match = 0;
        foreach $c (@common) {
          $p = ${$c}[0];
          $a = ${$c}[1];
          if ($keyw =~ /^$p$/) {
            if (! &array_lookup($a, $keyp)) {
              push @{$a}, $keyp;
            }
            $match = 1;
            last;
          }
        }
        if (! $match) {
          foreach $l (@list) {
            if (${$l}[0] eq $keyp) {
              $a = ${$l}[1];
              $match = 1;
              last;
            }
          }
          if (! $match) {
            push @list, [ $keyp, [@comments]];
          }
          else {
            foreach $text (@comments) {
              push @{$a}, $text;
            }
          }
        }
      }
      @comments = ();
    }
  }
  else {
    $text = $_;
    if (/(.*)\s\Q$mode\E\}\s*$/) {
      $text = $1;
      $mode = "";
    }
    $text =~ s/^\s{3}//;
    $text =~ s/\s+$//;
    push @comments, $text;
  }
}

if ($mode) {
  die "ERROR: still scanning for closing \" $mode}\"\n";
}


##############################################################################
# Print parameter names matched in common.idb.
##############################################################################

print "#" x 79, "\n";
print "# Parameter names matched in common.idb:\n";

foreach $c (@common) {
  foreach $keyp (@{${$c}[1]}) {
    print "#     $keyp\n";
  }
}

print "#" x 79, "\n";
print "#\n";


##############################################################################
# Print new keyword patterns.
##############################################################################

foreach $l (@list) {
  print "vars: ${$l}[0]\n";
  print "info:";
  $newl = 1;
  foreach $text (@{${$l}[1]}) {
    if ($newl) {
      print " ";
      $newl = 0;
    }
    print "$text\n";
  }
  print "\n" if ($newl);
  print "#\n";
}

__END__

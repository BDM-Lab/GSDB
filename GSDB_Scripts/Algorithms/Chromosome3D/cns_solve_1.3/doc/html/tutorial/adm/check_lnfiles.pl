#! /usr/bin/perl -w

@files = ();

use File::Find;

sub wanted {
  if (/^tutorial.csh$/) {
    push @files, $File::Find::name;
  }
}

find(\&wanted, '.');

foreach $f (@files) {
  print "$f\:\n";

  open INP, "< $f" or die "Cannot read file \"$f\"\n";

  $mode = 0;

  while (<INP>) {
    $file = "";
    if ($mode == 0) {
      if (/set\s+lnfiles\s*=\s*\((\s*[^\s\\\)]*)\s*(.?)/) {
        $file = $1;
        $mode = 1 if ($2 ne ")");
      }
    }
    else {
      if (/\s*(\S+)\s*([\\\)])\s*$/) {
        $file = $1;
        $mode = 0 if ($2 eq ")");
      }
      else {
        die "ERROR: improper line\n";
      }
    }

    if ($file) {
      if (! -f $file) {
        print "ERROR: $file\n";
      }
      else {
        print "   OK: $file\n";
      }
    }
  }

  close INP;
}

__END__

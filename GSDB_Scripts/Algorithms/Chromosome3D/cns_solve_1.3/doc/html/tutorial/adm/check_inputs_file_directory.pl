#! /usr/bin/perl -w

@files = ();

use File::Find;

sub wanted {
  return if ($File::Find::dir eq "./xtal_data");
  if (/^.*\.inp$/ or /^.*\.def$/) {
    push @files, $File::Find::name;
  }
}

find(\&wanted, '.');

foreach $f (@files) {
  open INP, "< $f" or die "Cannot read file \"$f\"\n";

  $master_directory = "";
  $master_file = "";

  while (<INP>) {
    if    (/\{\+ file:\s+(\S+)\s+\+\}/) {
      if ($master_file) {
        die "Duplicate {+ file: ... +} in \"$f\"\n";
      }
      $master_file = $1;
    }
    elsif (/\{\+ directory:\s+(\S+)\s+\+\}/) {
      if ($master_directory) {
        die "Duplicate {+ directory: ... +} in \"$f\"\n";
      }
      $master_directory = $1;
    }
  }

  close INP;

  if (not $master_file) {
    print "Missing {+ file: ... +} in \"$f\"\n";
  }
  if (not $master_directory) {
    print "Missing {+ directory: ... +} in \"$f\"\n";
  }
  if (    $master_file
      and $master_directory
      and $f ne "./$master_directory/$master_file") {
    print "Directory/file name mismatch\n";
    print "$f\n";
    print "./$master_directory/$master_file\n";
  }
}

__END__

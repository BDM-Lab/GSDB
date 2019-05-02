#! /usr/bin/perl -w

die "Error: environment variable CNS_SOLVE not defined.\n"
  if (not exists $ENV{CNS_SOLVE});

$CNS_INPUTS = "$ENV{CNS_SOLVE}/inputs";
$cns_transfer = "$ENV{CNS_SOLVE}/bin/cns_transfer";

@files = ();

use English;
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

  $master_path = "";

  if ($master_file and $master_directory) {
    $master_path = "$CNS_INPUTS/$master_directory/$master_file";
    if (not -f $master_path) {
      print "$f\:\n";
      print "Master file \"$master_path\" not found.\n";
    }
  }

  if ($master_path) {
    if (system(
      "$cns_transfer -def $f -inp $master_path > $f.trf") != 0) {
      die "cns_transfer failed.\n";
    }
    $diff_result = `diff -b -w -q $f $f.trf`;
    if ($CHILD_ERROR == 2) {
      die "diff failed.\n";
    }
    if ($diff_result ne "") {
      print "Changes: $f\n";
      system("twdiff $f $f.trf");
#      if (system(
#        "xdiff -geom 1100x850+150+100 -b -w $f $f.trf") != 0) {
#        die "xdiff failed.\n";
#      }
      $answer = `xmessage -geom +130+40 -buttons NO,yes 'Replace?' -print`;
      if ($CHILD_ERROR == 1) {
        die "xmessage failed.\n";
      }
      if ($answer =~ /yes\s*/) {
        if (system("mv $f $f.old") != 0) {
          die "ERROR: mv $f $f.old\n";
        }
        if (system("mv $f.trf $f") != 0) {
          die "ERROR: mv $f.trf $f\n";
        }
      }
    }
    else {
      if (system("rm $f.trf") != 0) {
        die "ERROR: rm $f.trf\n";
      }
      print "     OK: $f\n";
    }
  }
}

__END__

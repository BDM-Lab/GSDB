#!/usr/bin/perl
##############################################################################
#### Script to run CNS-DEN with various values of given parameters
#### and run them on a Linux cluster (using LFS program in this case).
#### This script can also be used to submit non-cluster jobs by changing a
#### parameter below.
####
#### Abhinav Kumar, Debanu Das
#### JCSG (JOINT CENTER FOR STRUCTURAL GENOMICS) www.jcsg.org
#### Oct 2011
##############################################################################
use warnings;
use strict;
use Cwd;
##############################################################################
#### Edit seven variables below
##############################################################################
my $use_cluster      = 'No';       ## change it to NO if not using cluster
my $dir_to_run_cns = "/Users/brunger/temp";    ## place to run cns
#### For a cluster submission, use scratch area
my $cv_file      = "/Users/brunger/Desktop/DEN/hp3342_data_n.cns";      ## absolute path for CNS reflection file
my $pdb_file     = "/Users/brunger/Desktop/DEN/phaser_shifted.pdb";    ## absolute path for coordinates (starting=reference coordinates)
my $ref_inp_file = "/Users/brunger/Desktop/DEN/refine_den.inp";           ## absolute path for refine_den.inp file
##### ref_inp_file should contain correct pdb and reflection file names
my %variables = (                              ### edit the list below, as well as add/subtract lines
#### 'parameters'  =>  [ 'val1', 'val2, ...],
   'den_gamma'      => [ '0.0', '0.2', '0.4', '0.6', '0.8', '1.0' ],
   'den_scale'      => [ '0', '1', '10', '30', '100', '300' ],
###   'minimize_nstep' => [ '10',   '100' ],
###   'temperature'    => [ '1000', '3000', '5000' ],
###   'seed'           => [ '1',    '5', '9' ],
);
my $cns_solve_setup = "/Users/Shared/cns_development/cns_solve_env";
##############################################################################
#### No edit below this line
##############################################################################
my $dir  = getcwd();
my $user = getlogin();
my @ref_inp_lines;
my $inp_dir = "$dir/Inputs";     ## store all inp files here
my $out_dir = "$dir/Outputs";    ## store all cns refined pdb and log files here
my $que_dir = "$dir/Que";        ## store all que log files for cluster here
##############################################################################
Make_dirs();                     ## make some directories
Read_inp($ref_inp_file);         ## read in the reference inp file
my @combinations = Generate_combinations();    ## generate all possible combinations
##############################################################################
foreach my $combination (@combinations) {
   print "Working on combination: $combination\n";
   my $new_inp  = "$inp_dir/refine_den_$combination.inp";
   my $cns_file = "$inp_dir/cns_$combination.sh";
   my $que_file = "$out_dir/que_$combination.log";
   Edit_file( $new_inp, $combination );
   Create_cns_script( $cns_file, $new_inp );
   Submit_jobs( $cns_file, $que_file );
}
print "Input script files are located in '$inp_dir'\n";
print "Output files and cluster job logs are located in '$out_dir'\n";
print "\n Keep checking your jobs with 'bjobs' \n";
print "\n Run 'sort_rfree.pl' after all jobs finish.\n";
##############################################################################
sub Create_cns_script {
   my ( $file, $inp ) = @_;
   my $id = $inp;
   $id =~ s/^.*refine_den_//;
   $id =~ s/\.inp//;
   open OUT, ">$file" || die "can't create file $file\n";
   print OUT "#!/bin/tcsh -f\n";
   print OUT "source $cns_solve_setup\n";
   print OUT "mkdir -p $dir_to_run_cns/$user/$id\n";
   print OUT "cd $dir_to_run_cns/$user/$id\n";
   print OUT "#rm -rf ./* \n";
   print OUT "cp $pdb_file $dir_to_run_cns/$user/$id/ \n";
   print OUT "cp $cv_file  $dir_to_run_cns/$user/$id/ \n";
   print OUT "cp $inp $dir_to_run_cns/$user/$id/run.inp \n";
   print OUT "cns_solve < run.inp > run.log \n";
   print OUT "cp cns_refined_new.pdb  $out_dir/refine_den_$id.pdb \n";
   print OUT "cp run.log $out_dir/refine_den_$id.log \n";
   close OUT;
   chmod 0755, $file;
}
##############################################################################
sub Edit_file {
   my ( $file, $combination ) = @_;
   open INP, ">$file" || die "can't create file $file\n";
   foreach my $value ( split /,/, $combination ) {
      my ( $param, $val ) = split /:/, $value;
      for (@ref_inp_lines) {
         s/^{===>} $param=.*$/{===>} $param=$val;/;
         s/^{===>} output_root=\".*$/{===>} output_root=\"cns_refined_new\";/;
      }
   }
### change the output root to 'cns_refined_new'
   for (@ref_inp_lines) {
      s/^{===>} output_root=\".*$/{===>} output_root=\"cns_refined_new\";/;
      print INP $_;
   }
   close INP;
}
##############################################################################
sub Generate_combinations {
   my @combinations;
   foreach my $item ( sort keys %variables ) {
      if (@combinations) {
         my @tmp;
         foreach my $entry (@combinations) {
            for ( my $i = 0 ; $i < @{ $variables{$item} } ; $i++ ) {
               push @tmp, "$entry,$item:$variables{$item}[$i]";
            }
         }
         @combinations = @tmp;
      }
      else {
         for ( my $i = 0 ; $i < @{ $variables{$item} } ; $i++ ) {
            push @combinations, "$item:$variables{$item}[$i]";
         }
      }
   }
   return @combinations;
}
##############################################################################
sub Make_dirs {
   for ( $inp_dir, $out_dir, $que_dir ) {
      mkdir $_ unless -e $_;
   }
}
##############################################################################
sub Read_inp {
   my $file = shift;
   open IN, $file || die "can't read file '$file'\n";
   @ref_inp_lines = <IN>;
   close IN;
}
##############################################################################
sub Submit_jobs {
   my ( $cns_file, $que_file ) = @_;
   if ( $use_cluster =~ /Yes/i ) {
      ### lsf que using que name sdcq; modify options as required
      system("bsub -q sdcq -o $que_file nice $cns_file");
   }
   else {
      system("nice $cns_file < /dev/null > cns_run.log 2>& 1 &");
   }
}
##############################################################################


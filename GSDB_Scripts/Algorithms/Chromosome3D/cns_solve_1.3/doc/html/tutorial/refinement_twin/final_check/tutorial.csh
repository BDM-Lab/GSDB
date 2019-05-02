#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/por/water_pick_twin.pdb \
	     data/por/porin.cv \
	     refinement_twin/final_check/omac \
             refinement_twin/final_check/model_stats_twin_final.inp \
	     refinement_twin/final_check/gradient_map_twin_final.inp \
	     refinement_twin/final_check/model_map_twin_final.inp)
#
if ( ! $?CNS_SOLVE ) then
  echo "CNS_SOLVE not defined"
  echo "CNS must be correctly setup to running this tutorial script"
  exit 1
endif
#
if ( ! -d $CNS_SOLVE/$tutdir ) then
  echo "cannot find tutorial directory $CNS_SOLVE/$tutdir"
  echo "please check that you are using the correct version of CNS"
  exit 1
endif
#
foreach file ( $lnfiles )
  if ( ! -e {$file}:t ) then
    echo "linking $tutdir/$file to current directory"
    ln -s $CNS_SOLVE/$tutdir/$file .
  endif
end
#
echo "        "
echo "running:"
echo "         cns_solve < model_stats_twin_final.inp > model_stats_twin_final.out"
echo "        "
#
cns_solve < model_stats_twin_final.inp > model_stats_twin_final.out
#
echo "        "
echo "running:"
echo "         cns_solve < gradient_map_twin_final.inp > gradient_map_twin_final.out"
echo "        "
#
cns_solve < gradient_map_twin_final.inp > gradient_map_twin_final.out
#
echo "        "
echo "running:"
echo "         cns_solve < model_map_twin_final.inp > model_map_twin_final.out"
echo "        "
#
cns_solve < model_map_twin_final.inp > model_map_twin_final.out
#
echo "        "
#

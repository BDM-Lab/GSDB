#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Axel Brunger
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(fret/fret_distribution_calculations/410_554.pdb \
             fret/fret_distribution_calculations/pseudo_atom_patches.def \
             fret/fret_distribution_calculations/model_anneal.inp \
             fret/fret_distribution_calculations/distance_calculation.inp \
             fret/fret_distribution_calculations/distance.def \
             fret/fret_distribution_calculations/model_anneal_fixed.inp )
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
echo "         cns_solve < model_anneal.inp > model_anneal.out"
echo "        "
#
cns_solve < model_anneal.inp > model_anneal.out
#
echo "        "
#
echo "running:"
echo "         cns_solve < distance_calculation.inp > distance_calculation.out"
echo "        "
#
cns_solve < distance_calculation.inp > distance_calculation.out
#
echo "        "
echo "running:"
echo "         cns_solve < model_anneal_fixed.inp > model_anneal_fixed.out"
echo "        "
#
cns_solve < model_anneal_fixed.inp > model_anneal_fixed.out

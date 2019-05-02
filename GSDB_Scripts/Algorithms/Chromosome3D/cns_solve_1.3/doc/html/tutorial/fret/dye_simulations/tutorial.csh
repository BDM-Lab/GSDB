#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Axel Brunger
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(fret/dye_simulations/140_350.pdb \
             fret/dye_simulations/cy3.pdb \
             fret/dye_simulations/cy5.pdb \
             fret/dye_simulations/dye_patches.def \
             fret/dye_simulations/model_anneal.inp \
             fret/dye_simulations/distance_calculation.inp \
             fret/dye_simulations/average_coordinates.inp)
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
echo "         cns_solve < average_coodinates.inp > average_coordinates.out"
echo "        "
#
cns_solve < average_coordinates.inp > average_coordinates.out
#

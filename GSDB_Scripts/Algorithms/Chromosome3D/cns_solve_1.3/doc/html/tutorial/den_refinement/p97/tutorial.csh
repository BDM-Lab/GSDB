#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/p97/p97_adp_start.pdb \
             data/p97/p97_adp_reference.pdb \
             data/p97/p97_adp_cns.par \
             data/p97/p97_adp_cns.top \
             data/p97/p97_adp_psad38.hkl \
             data/p97/p97_adp_testset_to38.hkl \
             data/p97/p97_adp_ncs.def \
             den_refinement/p97/refine_den.inp)
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
echo "         cns_solve < refine_den.inp > refine_den.out"
echo "        "
#
cns_solve < refine_den.inp > refine_den.out
#
echo "        "
#

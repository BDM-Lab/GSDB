#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/pdz/pdz.mtf \
             data/pdz/pdz.pdb \
             coord_manip/ncs_trans/merge_structures.inp \
             coord_manip/ncs_trans/realspace_transform.inp)
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
echo "         cns_solve < merge_structures.inp > merge_structures.out"
echo "        "
#
cns_solve < merge_structures.inp > merge_structures.out
#
echo "        "
echo "running:"
echo "         cns_solve < realspace_transform.inp > realspace_transform.out"
echo "        "
#
cns_solve < realspace_transform.inp > realspace_transform.out
#
echo "        "
#

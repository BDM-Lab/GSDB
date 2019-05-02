#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(generate/capping/generate.inp \
             generate/capping/capping.top \
             generate/capping/capping.param \
             generate/capping/protein_capping.link \
             generate/capping/capping.pdb \
             generate/capping/incomplete.inp \
             generate/capping/incomplete.param \
             generate/capping/incomplete.pdb)
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
echo "         cns_solve < generate.inp > generate.out"
echo "        "
#
               cns_solve < generate.inp > generate.out
#
echo "        "
echo "running:"
echo "         cns_solve < incomplete.inp > incomplete.out"
echo "        "
#
               cns_solve < incomplete.inp > incomplete.out
#
echo "        "
#

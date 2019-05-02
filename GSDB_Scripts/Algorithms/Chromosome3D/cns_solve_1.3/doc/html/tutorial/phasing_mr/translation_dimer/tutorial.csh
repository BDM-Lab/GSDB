#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/mbp/mbp_dimer_search.pdb \
             data/mbp/mbp_cns.pdb \
	     data/mbp/mbp_mr.hkl \
	     data/mbp/cross_rotation.list \
	     data/mbp/translation.pdb \
             phasing_mr/translation_dimer/translation_dimer.inp \
	     phasing_mr/translation_dimer/shift_molecules.inp)
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
#
echo "        "
echo "running:"
echo "         cns_solve < translation_dimer.inp > translation_dimer.out"
echo "        "
#
cns_solve < translation_dimer.inp > translation_dimer.out
#
echo "        "
echo "running:"
echo "         cns_solve < shift_molecules.inp > shift_molecules.out"
echo "        "
#
cns_solve < shift_molecules.inp > shift_molecules.out
#
echo "        "
#

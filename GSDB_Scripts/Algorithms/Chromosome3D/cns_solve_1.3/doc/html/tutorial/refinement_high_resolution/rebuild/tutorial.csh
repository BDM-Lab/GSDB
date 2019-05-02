#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/mbp/refine_yb.pdb \
	     data/mbp/mbp.cv \
	     data/mbp/mbp_o.pdb \
	     refinement_high_resolution/rebuild/sa_omit_map.inp \
	     refinement_high_resolution/rebuild/refine_reb.inp)
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
echo "         cns_solve < sa_omit_map.inp > sa_omit_map.out"
echo "        "
#
cns_solve < sa_omit_map.inp > sa_omit_map.out
#
echo "        "
echo "running:"
echo "         cns_solve < refine_reb.inp > refine_reb.out"
echo "        "
#
cns_solve < refine_reb.inp > refine_reb.out
#

#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/nsf/refine_atp.pdb \
	     data/nsf/nsf.cv \
	     data/nsf/mad_phase2.hkl \
	     data/nsf/atp.param \
	     data/nsf/atp.top \
	     refinement_mlhl/model_complete/refine_atp_2.inp )
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
echo "         cns_solve < refine_atp_2.inp > refine_atp_2.out"
echo "        "
#
cns_solve < refine_atp_2.inp > refine_atp_2.out
#
echo "        "
#

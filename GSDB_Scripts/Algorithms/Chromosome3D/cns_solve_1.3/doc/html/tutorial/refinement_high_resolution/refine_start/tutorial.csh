#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/mbp/mbp.pdb \
	     data/mbp/mbp.cv \
	     data/mbp/ncs_restrain.def \
             refinement_high_resolution/refine_start/rigid.inp \
             refinement_high_resolution/refine_start/model_map.inp \
             refinement_high_resolution/refine_start/refine.csh \
             refinement_high_resolution/refine_start/refine_ncs.inp \
	     refinement_high_resolution/refine_start/refine.inp)
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
echo "         cns_solve < rigid.inp > rigid.out"
echo "        "
#
cns_solve < rigid.inp > rigid.out
#
echo "        "
echo "running:"
echo "         cns_solve < refine.inp > refine.out"
echo "        "
#
cns_solve < refine.inp > refine.out
#
echo "        "
echo "running:"
echo "         cns_solve < refine_ncs.inp > refine_ncs.out"
echo "        "
#
cns_solve < refine_ncs.inp > refine_ncs.out
#

echo "        "
echo "running:"
echo "         ./refine.csh 1"
echo "        "
#
./refine.csh 1
#
echo "        "
echo "running:"
echo "         ./refine.csh 2"
echo "        "
#
./refine.csh 2
echo "        "
echo "running:"
echo "         cns_solve < model_map.inp > model_map.out"
echo "        "
#
cns_solve < model_map.inp > model_map.out

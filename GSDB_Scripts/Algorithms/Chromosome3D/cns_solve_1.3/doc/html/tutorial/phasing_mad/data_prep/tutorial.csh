#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/nsf/w1.sca \
             data/nsf/w2.sca \
             data/nsf/w3.sca \
             data/nsf/w4.sca \
             phasing_mad/data_prep/merge.inp \
             phasing_mad/data_prep/scale.inp)
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
echo "         to_cns w1.sca w1.hkl"
echo "        "
#
to_cns w1.sca w1.hkl
#
echo "        "
echo "running:"
echo "         to_cns w2.sca w2.hkl"
echo "        "
#
to_cns w2.sca w2.hkl
#
echo "        "
echo "running:"
echo "         to_cns w3.sca w3.hkl"
echo "        "
#
to_cns w3.sca w3.hkl
#
echo "        "
echo "running:"
echo "         to_cns w4.sca w4.hkl"
echo "        "
#
to_cns w4.sca w4.hkl
#
echo "        "
echo "running:"
echo "         cns_solve < merge.inp > merge.out"
echo "        "
#
cns_solve < merge.inp > merge.out
#
echo "        "
echo "running:"
echo "         cns_solve < scale.inp > scale.out"
echo "        "
#
cns_solve < scale.inp > scale.out
#
echo "        "
#

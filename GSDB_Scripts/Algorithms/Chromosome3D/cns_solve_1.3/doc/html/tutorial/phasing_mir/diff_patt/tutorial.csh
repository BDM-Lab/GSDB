#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/pen/scale.hkl \
             data/pen/kuof_sites.sdb \
             phasing_mir/diff_patt/patterson_map.inp \
             phasing_mir/diff_patt/predict_patterson.inp)
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
echo "         cns_solve < patterson_map.inp > patterson_map.out"
echo "        "
#
cns_solve < patterson_map.inp > patterson_map.out
#
echo "        "
echo "running:"
echo "         plot_patterson"
echo "        "
#
plot_patterson << EOF
patterson_map_y.map
y
patterson_map_y.ps



Experimental Patterson map Y-section
EOF
#
echo "        "
echo "running:"
echo "         cns_solve < predict_patterson.inp > predict_patterson.out"
echo "        "
#
cns_solve < predict_patterson.inp > predict_patterson.out
#
echo "        "
echo "running:"
echo "         plot_patterson"
echo "        "
#
plot_patterson << EOF
predict_patterson_y.map
y
predict_patterson_y.ps



Predicted Patterson map Y-section
EOF
#
echo "        "
#

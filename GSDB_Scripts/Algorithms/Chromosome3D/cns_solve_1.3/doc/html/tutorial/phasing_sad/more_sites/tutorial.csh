#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams & Ralf W. Grosse-Kunstleve
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/nsf/sad_scale.hkl \
             data/nsf/sad_phase.hkl \
             data/nsf/sad_phase_grad.hkl \
             data/nsf/sad_phase.sdb \
             phasing_sad/more_sites/fourier_map_anom.inp \
             phasing_sad/more_sites/fourier_map_grad.inp \
             phasing_sad/more_sites/sdb_manipulate.inp \
             phasing_sad/more_sites/omac)
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
echo "         cns_solve < fourier_map_anom.inp > fourier_map_anom.out"
echo "        "
#
               cns_solve < fourier_map_anom.inp > fourier_map_anom.out
#
echo "        "
echo "running:"
echo "         cns_solve < fourier_map_grad.inp > fourier_map_grad.out"
echo "        "
#
               cns_solve < fourier_map_grad.inp > fourier_map_grad.out
#
echo "        "
echo "running:"
echo "         cns_solve < sdb_manipulate.inp > sdb_manipulate.out"
echo "        "
#
               cns_solve < sdb_manipulate.inp > sdb_manipulate.out
#
echo "        "
#

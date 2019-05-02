#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams & Ralf W. Grosse-Kunstleve
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/nsf/mad_scale.hkl \
             data/nsf/mad_phase_flip.hkl \
             data/nsf/mad_phase_flip_grad.hkl \
             data/nsf/mad_phase_flip.sdb \
             data/nsf/mad_eight_sites.sdb \
             phasing_mad/more_sites/fourier_map_anom_flip.inp \
             phasing_mad/more_sites/fourier_map_grad_flip.inp \
             phasing_mad/more_sites/sdb_manipulate.inp \
             phasing_mad/more_sites/mad_phase2.inp \
             phasing_mad/more_sites/omac)
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
echo "         cns_solve < fourier_map_anom_flip.inp > fourier_map_anom_flip.out"
echo "        "
#
               cns_solve < fourier_map_anom_flip.inp > fourier_map_anom_flip.out
#
echo "        "
echo "running:"
echo "         cns_solve < fourier_map_grad_flip.inp > fourier_map_grad_flip.out"
echo "        "
#
               cns_solve < fourier_map_grad_flip.inp > fourier_map_grad_flip.out
#
echo "        "
echo "running:"
echo "         cns_solve < sdb_manipulate.inp > sdb_manipulate.out"
echo "        "
#
               cns_solve < sdb_manipulate.inp > sdb_manipulate.out
#
echo "        "
echo "running:"
echo "         cns_solve < mad_phase2.inp > mad_phase2.out"
echo "        "
#
               cns_solve < mad_phase2.inp > mad_phase2.out
#
echo "        "
#

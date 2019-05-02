#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams & Ralf W. Grosse-Kunstleve
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/nsf/mad_scale.hkl \
             data/nsf/heavy_search_1.sdb \
             phasing_mad/heavy_refine/flip_sites.inp \
             phasing_mad/heavy_refine/mad_phase.inp \
             phasing_mad/heavy_refine/mad_phase_flip.inp)
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
echo "         cns_solve < flip_sites.inp > flip_sites.out"
echo "        "
#
               cns_solve < flip_sites.inp > flip_sites.out
#
echo "        "
echo "running:"
echo "         cns_solve < mad_phase.inp > mad_phase.out"
echo "        "
#
               cns_solve < mad_phase.inp > mad_phase.out
#
echo "        "
echo "running:"
echo "         cns_solve < mad_phase_flip.inp > mad_phase_flip.out"
echo "        "
#
               cns_solve < mad_phase_flip.inp > mad_phase_flip.out
#
echo "        "
#

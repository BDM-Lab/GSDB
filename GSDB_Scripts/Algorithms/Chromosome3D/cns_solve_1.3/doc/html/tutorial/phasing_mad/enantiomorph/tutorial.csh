#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams & Ralf W. Grosse-Kunstleve
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/nsf/mad_scale.hkl \
             data/nsf/mad_phase.hkl \
             data/nsf/mad_phase_flip.hkl \
             phasing_mad/enantiomorph/fourier_map_elec.inp \
             phasing_mad/enantiomorph/fourier_map_elec_flip.inp \
             phasing_mad/enantiomorph/omac)
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
echo "         cns_solve < fourier_map_elec.inp > fourier_map_elec.out"
echo "        "
#
               cns_solve < fourier_map_elec.inp > fourier_map_elec.out
#
echo "        "
echo "running:"
echo "         cns_solve < fourier_map_elec_flip.inp > fourier_map_elec_flip.out"
echo "        "
#
               cns_solve < fourier_map_elec_flip.inp > fourier_map_elec_flip.out
#
echo "        "
#

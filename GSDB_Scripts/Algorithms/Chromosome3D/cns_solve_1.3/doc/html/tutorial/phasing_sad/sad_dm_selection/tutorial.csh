#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams & Ralf W. Grosse-Kunstleve
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/nsf/sad_scale.hkl \
             data/nsf/sad_eight_sites.sdb \
             phasing_sad/sad_dm_selection/flip_sites.inp \
             phasing_sad/sad_dm_selection/sad_phase2.inp \
             phasing_sad/sad_dm_selection/sad_phase2_flip.inp \
             phasing_sad/sad_dm_selection/matthews_coef.inp \
             phasing_sad/sad_dm_selection/density_modify.inp \
             phasing_sad/sad_dm_selection/density_modify_flip.inp \
             phasing_sad/sad_dm_selection/fourier_map_sad2.inp \
             phasing_sad/sad_dm_selection/fourier_map_sad2_flip.inp \
             phasing_sad/sad_dm_selection/omac)
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
echo "         cns_solve < sad_phase2.inp > sad_phase2.out"
echo "        "
#
               cns_solve < sad_phase2.inp > sad_phase2.out
#
echo "        "
echo "running:"
echo "         cns_solve < sad_phase2_flip.inp > sad_phase2_flip.out"
echo "        "
#
               cns_solve < sad_phase2_flip.inp > sad_phase2_flip.out
#
echo "        "
echo "running:"
echo "         cns_solve < matthews_coef.inp > matthews_coef.out"
echo "        "
#
               cns_solve < matthews_coef.inp > matthews_coef.out
#
echo "        "
echo "running:"
echo "         cns_solve < density_modify.inp > density_modify.out"
echo "        "
#
               cns_solve < density_modify.inp > density_modify.out
#
echo "        "
echo "running:"
echo "         cns_solve < density_modify_flip.inp > density_modify_flip.out"
echo "        "
#
               cns_solve < density_modify_flip.inp > density_modify_flip.out
#
echo "        "
echo "running:"
echo "         cns_solve < fourier_map_sad2.inp > fourier_map_sad2.out"
echo "        "
#
               cns_solve < fourier_map_sad2.inp > fourier_map_sad2.out
#
echo "        "
echo "running:"
echo "         cns_solve < fourier_map_sad2_flip.inp > fourier_map_sad2_flip.out"
echo "        "
#
               cns_solve < fourier_map_sad2_flip.inp > fourier_map_sad2_flip.out
#
echo "        "
#

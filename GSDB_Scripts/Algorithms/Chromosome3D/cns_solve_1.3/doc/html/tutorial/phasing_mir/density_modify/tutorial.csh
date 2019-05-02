#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams & Ralf W. Grosse-Kunstleve
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/pen/scale.hkl \
             data/pen/mir_phase_final.hkl \
             data/pen/mir_phase_final.sdb \
             phasing_mir/density_modify/matthews_coef.inp \
             phasing_mir/density_modify/density_modify.inp \
             phasing_mir/density_modify/fourier_map_mir.inp \
             phasing_mir/density_modify/fourier_map_dm.inp \
             phasing_mir/density_modify/sdb_to_pdb.inp \
             phasing_mir/density_modify/omac)
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
echo "         cns_solve < sdb_to_pdb.inp > sdb_to_pdb.out"
echo "        "
#
               cns_solve < sdb_to_pdb.inp > sdb_to_pdb.out
#
echo "        "
echo "running:"
echo "         cns_solve < fourier_map_mir.inp > fourier_map_mir.out"
echo "        "
#
               cns_solve < fourier_map_mir.inp > fourier_map_mir.out
#
echo "        "
echo "running:"
echo "         cns_solve < fourier_map_dm.inp > fourier_map_dm.out"
echo "        "
#
               cns_solve < fourier_map_dm.inp > fourier_map_dm.out
#
echo "        "
#

#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams & Ralf W. Grosse-Kunstleve
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/pen/scale.hkl \
             data/pen/ir_phase_kuof.hkl \
             data/pen/kuof_sites.sdb \
	     data/pen/phga.pdb \
             phasing_mir/cross_phase/cross_phase_phga.inp \
             phasing_mir/cross_phase/mir_phase_2deriv.inp \
             phasing_mir/cross_phase/pdb_to_sdb.inp \
             phasing_mir/cross_phase/omac)
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
echo "         cns_solve < cross_phase_phga.inp > cross_phase_phga.out"
echo "        "
#
               cns_solve < cross_phase_phga.inp > cross_phase_phga.out
#
echo "        "
echo "running:"
echo "         cns_solve < pdb_to_sdb.inp > pdb_to_sdb.out"
echo "        "
#
               cns_solve < pdb_to_sdb.inp > pdb_to_sdb.out
#
echo "        "
echo "running:"
echo "         cns_solve < mir_phase_2deriv.inp > mir_phase_2deriv.out"
echo "        "
#
               cns_solve < mir_phase_2deriv.inp > mir_phase_2deriv.out
#
echo "        "
#

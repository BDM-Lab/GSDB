#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams & Ralf W. Grosse-Kunstleve
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/pen/scale.hkl \
             data/pen/mir_phase_2deriv.hkl \
             data/pen/mir_phase_2deriv_grad.hkl \
             data/pen/kuof_sites.sdb \
             data/pen/phga_sites.sdb \
	     data/pen/kuof_sites_all.sdb \
	     data/pen/phga_sites_all.sdb \
             phasing_mir/more_sites/fourier_phga.inp \
	     phasing_mir/more_sites/gradient_phga.inp \
             phasing_mir/more_sites/sdb_manipulate.inp \
             phasing_mir/more_sites/mir_phase_final.inp \
             phasing_mir/more_sites/omac)
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
echo "         cns_solve < fourier_phga.inp > fourier_phga.out"
echo "        "
#
               cns_solve < fourier_phga.inp > fourier_phga.out
#
echo "        "
echo "running:"
echo "         cns_solve < gradient_phga.inp > gradient_phga.out"
echo "        "
#
               cns_solve < gradient_phga.inp > gradient_phga.out
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
echo "         cns_solve < mir_phase_final.inp > mir_phase_final.out"
echo "        "
#
               cns_solve < mir_phase_final.inp > mir_phase_final.out
#
echo "        "
#

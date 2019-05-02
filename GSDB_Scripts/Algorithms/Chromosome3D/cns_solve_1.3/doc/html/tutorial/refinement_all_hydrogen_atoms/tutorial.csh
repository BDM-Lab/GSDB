#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/mbp/mbp.pdb \
	     data/mbp/mbp.cv \
	     data/mbp/ncs_restrain.def \
             data/mbp/mbp_water.pdb \
             refinement_all_hydrogen_atoms/generate_easy.inp \
             refinement_all_hydrogen_atoms/cis_peptide.inp \
             refinement_all_hydrogen_atoms/refine.inp )
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
echo "         cns_solve < generate_easy.inp > generate_easy.out"
echo "        "
#
cns_solve < generate_easy.inp > generate_easy.out
#
echo "        "
echo "running:"
echo "         cns_solve < cis_peptide.inp > cis_peptide.out"
echo "        "
#
cns_solve < cis_peptide.inp > cis_peptide.out
echo "        "
echo "        "
echo "running:"
echo "         cns_solve < refine.inp > refine.out"
echo "        "
#
cns_solve < refine.inp > refine.out
#

#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Paul Adams
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(data/pdb/2QAL.pdb \
             data/pdb/2QAM.pdb \
             data/pdb/2QAN.pdb \
             data/pdb/2QAO.pdb \
             data/pdb/2qao-sf.cif \
             data/pdb/2qao-sf.cif.CNS \
             generate/ribosome/generate_easy_2QAL.inp \
             generate/ribosome/generate_easy_2QAM.inp \
             generate/ribosome/generate_easy_2QAN.inp \
             generate/ribosome/generate_easy_2QAO.inp \
             generate/ribosome/nmy.top \
             generate/ribosome/nmy.param \
             generate/ribosome/merge.inp \
             generate/ribosome/refine.inp )
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
echo "         cns_solve < generate_easy_2QAL.inp > generate_easy_2QAL.out "
echo "        "
#
cns_solve < generate_easy_2QAL.inp > generate_easy_2QAL.out
#
echo "        "
#
echo "running:"
echo "         cns_solve < generate_easy_2QAM.inp > generate_easy_2QAM.out "
echo "        "
#
cns_solve < generate_easy_2QAM.inp > generate_easy_2QAM.out
#
echo "        "
echo "running:"
echo "         cns_solve < generate_easy_2QAN.inp > generate_easy_2QAN.out "
echo "        "
#
cns_solve < generate_easy_2QAN.inp > generate_easy_2QAN.out
#
echo "        "
echo "running:"
echo "         cns_solve < generate_easy_2QAO.inp > generate_easy_2QAO.out "
echo "        "
#
cns_solve < generate_easy_2QAO.inp > generate_easy_2QAO.out
#
echo "        "
echo "running:"
echo "         cns_solve < merge.inp > merge.out "
echo "        "
#
cns_solve < merge.inp > merge.out
#
echo "        "
echo "running:"
echo "         cns_solve < refine.inp > refine.out "
echo "        "
#
cns_solve < refine.inp > refine.out
#
echo "        "

#!/bin/csh -f
#
# tutorial script for CNSsolve
#
# written by: Axel Brunger
# copyright Yale University
#
set tutdir='doc/html/tutorial'
set lnfiles=(fret/docking_calculations/model_A_S61C.pdb \
             fret/docking_calculations/generate_easy_A_S61C_cy3_built.pdb \
             fret/docking_calculations/fret_start_model.pdb \
             fret/docking_calculations/pseudo_atom_patches.def \
             fret/docking_calculations/fret_restraints.def \
             fret/docking_calculations/generate_easy_A_S61C.inp \
             fret/docking_calculations/dye_patches.def \
             fret/docking_calculations/model_anneal_A_S61C.inp \
             fret/docking_calculations/average_coordinates_A_S61C.inp \
             fret/docking_calculations/docking.inp \
             fret/docking_calculations/rms_matrix.inp \
             fret/docking_calculations/sort.pl)
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
echo "         cns_solve < generate_easy_A_S61C.inp > generate_easy_A_S61C.out"
echo "        "
#
cns_solve < generate_easy_A_S61C.inp > generate_easy_A_S61C.out
#
echo "        "
#
echo "running:"
echo "         cns_solve < model_anneal_A_S61C.inp > model_anneal_A_S61C.out"
echo "        "
#
cns_solve < model_anneal_A_S61C.inp > model_anneal_A_S61C.out
#
echo "        "
echo "running:"
echo "         cns_solve < average_coodinates_A_S61C.inp > average_coordinates_A_S61C.out"
echo "        "
#
cns_solve < average_coordinates_A_S61C.inp > average_coordinates_A_S61C.out
#
echo "        "
echo "running:"
echo "         cns_solve < dockng.inp > docking.out"
echo "        "
#
cns_solve < docking.inp > docking.out
#
echo "        "
echo "running:"
echo "         ./sort.pl docking_*.pdb > sort.list"
echo "        "
#
./sort.pl docking_*.pdb > sort.list
#
echo "        "
echo "        "
echo "running:"
echo "         cns_solve < rms_matrix.inp > rms_matrix.out"
echo "        "
#
cns_solve < rms_matrix.inp > rms_matrix.out
#
echo "        "
echo "        "
echo "running:"
echo "         cluster_struc rms_matrix.list 3. 2. > cluster_struc.list"
echo "        "
#
cluster_struc rms_matrix.list 3. 2. > cluster_struc.list


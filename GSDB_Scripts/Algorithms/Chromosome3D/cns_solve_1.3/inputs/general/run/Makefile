###############################################################################
#									      #
#       Makefile for CNSsolve input run directory			      #
#       Paul Adams 2-9-98						      #
#       copyright Yale University                                             #
#									      #
###############################################################################

.SUFFIXES:
.SUFFIXES: .sed

SHELL=/bin/sh

FIRSTFILES=

FILES=alternate \
      bfactor_plot \
      cis_peptide \
      cns_to_o \
      delta_phipsi \
      rmsd \
      difference_distance \
      generate \
      generate_easy \
      generate_seq \
      get_ncs_matrices \
      merge_structures \
      model_mask \
      model_minimize \
      model_anneal \
      model_rigid \
      ramachandran \
      split_structure \
      surface_plot \
      buried_surface \
      contact \
      shift_molecules \
      realspace_transform \
      fractional_transform \
      neighbours \
      psf_to_mtf \
      rename_segid \
      delete_atoms \
      mtf_to_psf \
      molecule_extent \
      hydrogen_bonds \
      rms_fit \
      den_make \
      average_coordinates \
      distance_calculation \
      rms_matrix \
      docking

LINKS=../../xtal_data \
      ../../auxiliary \
      ..

include ../../Makefile.include

alternate.sed:

bfactor_plot.sed:

cns_to_o.sed:

cis_peptide.sed:

delta_phipsi.sed:

rmsd.sed:

difference_distance.sed:

generate.sed:

generate_easy.sed:
	@ sed -e \
	  "s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"2BF1.pdb\"/; \
           s/{===>} structure_outfile=\".*\"/{===>} structure_outfile=\"OUTPUT\"/; \
           s/{===>} atom_build=(none)/{===>} atom_build=(not known)/; \
           s/{===>} atom_delete=(not known)/{===>} atom_delete=(none)/; \
           s/{===>} coordinate_outfile=\".*\"/{===>} coordinate_outfile=\"generate_easy_1.pdb\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"3HXM.pdb\"/; \
           s/{===>} structure_outfile=\".*\"/{===>} structure_outfile=\"OUTPUT\"/; \
           s/{===>} coordinate_outfile=\".*\"/{===>} coordinate_outfile=\"generate_easy_2.pdb\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)
	@ sed -e \
	  "s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"3HXM_break.pdb\"/; \
           s/{===>} structure_outfile=\".*\"/{===>} structure_outfile=\"OUTPUT\"/; \
           s/{===>} coordinate_outfile=\".*\"/{===>} coordinate_outfile=\"generate_easy_3.pdb\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_3.inp)
	@ sed -e \
	  "s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"1AXY.pdb\"/; \
           s/{===>} structure_outfile=\".*\"/{===>} structure_outfile=\"OUTPUT\"/; \
           s/{===>} coordinate_outfile=\".*\"/{===>} coordinate_outfile=\"generate_easy_4.pdb\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_4.inp)
	@ sed -e \
	  "s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"2WN9.pdb\"/; \
           s/{===>} structure_outfile=\".*\"/{===>} structure_outfile=\"OUTPUT\"/; \
           s/{===>} coordinate_outfile=\".*\"/{===>} coordinate_outfile=\"generate_easy_5.pdb\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_5.inp)
	@ sed -e \
	  "s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"2WR0.pdb\"/; \
           s/{===>} structure_outfile=\".*\"/{===>} structure_outfile=\"OUTPUT\"/; \
           s/{===>} carbo_dist=.*/{===>} carbo_dist=3.2;/; \
           s/{===>} coordinate_outfile=\".*\"/{===>} coordinate_outfile=\"generate_easy_6.pdb\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_6.inp)
	@ sed -e \
	  "s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"2BF1_partial.pdb\"/; \
           s/{===>} structure_outfile=\".*\"/{===>} structure_outfile=\"OUTPUT\"/; \
           s/{===>} atom_build=(none)/{===>} atom_build=(not known)/; \
           s/{===>} atom_delete=(not known)/{===>} atom_delete=(none)/; \
           s/{===>} coordinate_outfile=\".*\"/{===>} coordinate_outfile=\"generate_easy_7.pdb\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_7.inp)



generate_seq.sed:

get_ncs_matrices.sed:

model_mask.sed:

model_minimize.sed:

model_anneal.sed:

model_rigid.sed:

merge_structures.sed:

ramachandran.sed:

split_structure.sed:

surface_plot.sed:

buried_surface.sed:

contact.sed:

shift_molecules.sed:

realspace_transform.sed:

fractional_transform.sed:

neighbours.sed:

psf_to_mtf.sed:

mtf_to_psf.sed:

rename_segid.sed:

delete_atoms.sed:

molecule_extent.sed:

hydrogen_bonds.sed:

rms_fit.sed:

den_make.sed:

average_coordinates.sed:

distance_calculation.sed:

rms_matrix.sed:

docking.sed:





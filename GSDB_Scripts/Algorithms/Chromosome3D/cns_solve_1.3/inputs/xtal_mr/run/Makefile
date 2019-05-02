###############################################################################
#									      #
#       Makefile for CNSsolve input run directory			      #
#       Paul Adams 2-9-98						      #
#	copyright Yale University					      #
#									      #
###############################################################################

.SUFFIXES:
.SUFFIXES: .sed

SHELL=/bin/sh

FIRSTFILES=cross_rotation

FILES=self_rotation \
      translation

LINKS=../../xtal_data \
      ../../auxiliary \
      ..

include ../../Makefile.include

cross_rotation.sed:
	@ sed -e \
	  "s/{===>} rot_method=\".*\"/{===>} rot_method=\"realspace\"/; \
	   s/{===>} rf_list_outfile=\".*\"/{===>} rf_list_outfile=\"cross_rotation_1.list\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

self_rotation.sed:

translation.sed:
	@ sed -e \
	  "s/{===>} reflection_infile_1=\".*\"/{===>} reflection_infile_1=\"3r9a_det.phase\"/; \
	   s/{===>} obs_phase=\".*\"/obs_phase=\"fobs\"/; \
	   s/{===>} obs_fom=\".*\"/obs_fom=\"fom\"/; \
	   s/{===>} pc_target=\".*\"/{===>} pc_target=\"f2f2\"/; \
	   s/{===>} rigid_target=\".*\"/{===>} rigid_target=\"vector\"/; \
	   s/{===>} trans_method=\".*\"/{===>} trans_method=\"phased\"/; \
	   s/{===>} ntrials=[0-9]*/{===>} ntrials=10/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"translation_phased\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"cross_trans.pdb\"/; \
	   s/{===>} atom_placed=(none)/{===>} atom_placed=( segid H )/; \
	   s/{===>} pcgroup_2=(resid 1:116 and segid H)/{===>} pcgroup_2=(resid 107:212 and segid L)/; \
	   s/{===>} pcgroup_3=(resid 107:212 and segid L)/{===>} pcgroup_3=(none)/; \
	   s/{===>} pcgroup_4=(resid 117:213 and segid H)/{===>} pcgroup_4=(none)/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"translation_cross\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)

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

FIRSTFILES=

FILES=make_cv_twin \
      twin_fraction \
      rigid_twin \
      bdomain_twin \
      anneal_twin \
      minimize_twin \
      bgroup_twin \
      bindividual_twin \
      model_map_twin \
      model_stats_twin \
      water_pick_twin \
      detwin_partial \
      detwin_perfect \
      detect_twinning \
      refine_twin \
      refine_twin_den \
      gradient_map_twin

LINKS=../../xtal_data \
      ../../auxiliary \
      ..

include ../../Makefile.include

make_cv_twin.sed:

twin_fraction.sed:

rigid_twin.sed:

bdomain_twin.sed:

anneal_twin.sed:
	@ sed -e \
	  "s/{===>} output_root=\".*\"/{===>} output_root=\"anneal_twin_den\"/; \
	   s/{===>} use_den=false/{===>} use_den=true/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

minimize_twin.sed:
	@ sed -e \
	  "s/{===>} coordinate_outfile=\".*\"/{===>} coordinate_outfile=\"minimize_twin_den.pdb\"/; \
	   s/{===>} use_den=false/{===>} use_den=true/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

bgroup_twin.sed:

bindividual_twin.sed:

model_map_twin.sed:

model_stats_twin.sed:

water_pick_twin.sed:

detwin_partial.sed:

detwin_perfect.sed:

detect_twinning.sed:

refine_twin.sed:
	@ sed -e \
	  "s/{===>} output_root=\".*\"/{===>} output_root=\"refine_twin_den\"/; \
           s/{===>} structure_infile=\".*\"/{===>} structure_infile=\"\"/; \
	   s/{===>} full_vdw_cycles=3/{===>} full_vdw_cycles=2/; \
	   s/{===>} den_infile=\".*\"/{===>} den_infile=\"porin.den\"/; \
	   s/{===>} use_den=false/{===>} use_den=true/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

refine_twin_den.sed:

gradient_map_twin.sed:

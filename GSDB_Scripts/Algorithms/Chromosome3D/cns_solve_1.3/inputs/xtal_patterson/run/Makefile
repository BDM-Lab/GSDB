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

FIRSTFILES=heavy_search

FILES=patterson_map \
      patterson_refine \
      predict_patterson

LINKS=../../xtal_data \
      ../../auxiliary \
      ..

include ../../Makefile.include

heavy_search.sed:

patterson_map.sed:
	@ sed -e \
	  "s/{===>} number_of_maps=[0-9]*/{===>} number_of_maps=3/; \
	   s/{===>} map.1.mode=\".*\"/{===>} map.1.mode=\"dispersive\"/; \
	   s/{===>} map.1.f_a=\".*\"/{===>} map.1.f_a=\"f_w4\"/; \
	   s/{===>} map.1.s_a=\".*\"/{===>} map.1.s_a=\"s_w4\"/; \
	   s/{===>} map.1.f_b=\".*\"/{===>} map.1.f_b=\"f_w1\"/; \
	   s/{===>} map.1.s_b=\".*\"/{===>} map.1.s_b=\"s_w1\"/; \
	   s/{===>} map.1.cut_df=0./{===>} map.1.cut_df=0.2/; \
	   s/{===>} map.1.max_df=100/{===>} map.1.max_df=6/; \
	   s/{===>} map.2.mode=\".*\"/{===>} map.2.mode=\"isomorphous\"/; \
	   s/{===>} map.2.f_a=\".*\"/{===>} map.2.f_a=\"f_w3\"/; \
	   s/{===>} map.2.s_a=\".*\"/{===>} map.2.s_a=\"s_w3\"/; \
	   s/{===>} map.2.f_b=\".*\"/{===>} map.2.f_b=\"f_w1\"/; \
	   s/{===>} map.2.s_b=\".*\"/{===>} map.2.s_b=\"s_w1\"/; \
	   s/{===>} map.3.mode=\".*\"/{===>} map.3.mode=\"anomalous\"/; \
	   s/{===>} map.3.f_a=\".*\"/{===>} map.3.f_a=\"f_w2\"/; \
	   s/{===>} map.3.s_a=\".*\"/{===>} map.3.s_a=\"s_w2\"/; \
           s/{===>} show_info=\".*\"/{===>} show_info=\"all\"/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"patterson_map_1\"/"  \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

patterson_refine.sed:

predict_patterson.sed:


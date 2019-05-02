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

FIRSTFILES=model_phase \
           merge \
           scale \
           analyse

FILES=combine \
      fourier_map \
      hlcoeff_blur \
      make_cv \
      make_hlcoeff \
      manipulate \
      matthews_coef \
      transform_map \
      average_friedels \
      flip_friedels \
      model_fcalc \
      average_map \
      mask_map \
      em_map_to_hkl \

LINKS=../../xtal_data \
      ../../auxiliary \
      ..

include ../../Makefile.include

analyse.sed:
	@ sed -e \
	  "s/{===>} ehistogram=false/{===>} ehistogram=true/; \
	   s/{===>} phases=false/{===>} phases=true/; \
	   s/{===>} fom=\".*\"/fom=\"fom\"/; \
	   s/{===>} reflection_infile_1=\".*\"/{===>} reflection_infile_1=\"model_phase.hkl\"/; \
	   s/{===>} reflection_infile_2=\".*\"/{===>} reflection_infile_2=\"amy_mir.hkl\"/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"analyse_1\"/; \
	   s/{===>} f_1=\".*\"/{===>} f_1=\"f_model\"/; \
	   s/{===>} f_2=\".*\"/{===>} f_2=\"fobs\"/; \
	   s/{===>} f_3=\".*\"/{===>} f_3=\"\"/; \
	   s/{===>} f_4=\".*\"/{===>} f_4=\"\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

combine.sed:

fourier_map.sed:
	@ sed -e \
	  "s/{===>} number_of_maps=[0-9]*/{===>} number_of_maps=3/; \
	   s/{===>} mode_1=\".*\"/{===>} mode_1=\"dispersive\"/; \
	   s/{===>} f_a_1=\".*\"/{===>} f_a_1=\"f_w4\"/; \
	   s/{===>} s_a_1=\".*\"/{===>} s_a_1=\"s_w4\"/; \
	   s/{===>} f_b_1=\".*\"/{===>} f_b_1=\"f_w1\"/; \
	   s/{===>} s_b_1=\".*\"/{===>} s_b_1=\"s_w1\"/; \
	   s/{===>} pa_1=\".*\"/{===>} pa_1=\"pa\"/; \
	   s/{===>} pb_1=\".*\"/{===>} pb_1=\"pb\"/; \
	   s/{===>} pc_1=\".*\"/{===>} pc_1=\"pc\"/; \
	   s/{===>} pd_1=\".*\"/{===>} pd_1=\"pd\"/; \
	   s/{===>} cut_df_1=0./{===>} cut_df_1=0.2/; \
	   s/{===>} max_df_1=100/{===>} max_df_1=6/; \
	   s/{===>} mode_2=\".*\"/{===>} mode_2=\"isomorphous\"/; \
	   s/{===>} f_a_2=\".*\"/{===>} f_a_2=\"f_w3\"/; \
	   s/{===>} s_a_2=\".*\"/{===>} s_a_2=\"s_w3\"/; \
	   s/{===>} f_b_2=\".*\"/{===>} f_b_2=\"f_w1\"/; \
	   s/{===>} s_b_2=\".*\"/{===>} s_b_2=\"s_w1\"/; \
	   s/{===>} pa_2=\".*\"/{===>} pa_2=\"pa\"/; \
	   s/{===>} pb_2=\".*\"/{===>} pb_2=\"pb\"/; \
	   s/{===>} pc_2=\".*\"/{===>} pc_2=\"pc\"/; \
	   s/{===>} pd_2=\".*\"/{===>} pd_2=\"pd\"/; \
	   s/{===>} mode_3=\".*\"/{===>} mode_3=\"anomalous\"/; \
	   s/{===>} f_a_3=\".*\"/{===>} f_a_3=\"f_w2\"/; \
	   s/{===>} s_a_3=\".*\"/{===>} s_a_3=\"s_w2\"/; \
	   s/{===>} pa_3=\".*\"/{===>} pa_3=\"\"/; \
	   s/{===>} pb_3=\".*\"/{===>} pb_3=\"\"/; \
	   s/{===>} pc_3=\".*\"/{===>} pc_3=\"\"/; \
	   s/{===>} pd_3=\".*\"/{===>} pd_3=\"\"/; \
	   s/{===>} phase_array_3=\".*\"/{===>} phase_array_3=\"mad\"/; \
	   s/{===>} fom_name_3=\".*\"/{===>} fom_name_3=\"fom\"/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"fourier_map_1\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/{===>} map_mode=\".*\"/{===>} map_mode=\"molecule\"/; \
	   s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"1YTT.pdb\"/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"fourier_map_2\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)

hlcoeff_blur.sed:

make_cv.sed:

make_hlcoeff.sed:

manipulate.sed:
	@ sed -e \
	  "s/new\.1\.name=\".*\"/new\.1\.name=\"new\"/; \
	   s/new\.2\.name=\".*\"/new\.2\.name=\"hlnew\"/; \
	   s/new\.2\.name=\".*\"/new\.2\.name=\"hlnew\"/; \
	   s/new\.2\.type=\".*\"/new\.2\.type=\"hendrickson-lattman\"/; \
	   s/new\.2\.type=\".*\"/new\.2\.type=\"hendrickson-lattman\"/; \
	   s/new\.1\.type=\".*\"/new\.1\.type=\"complex\"/; \
	   s/new\.1\.type=\".*\"/new\.1\.type=\"complex\"/; \
	   s/op\.1\.expression=.*;/op\.1\.expression=(new=fom*mad);/; \
	   s/op\.1\.selection=.*;/op\.1\.selection=(d>4);/; \
	   s/op\.2\.expression=.*;/op\.2\.expression=(hlnewa=pa);/; \
	   s/op\.2\.selection=.*;/op\.2\.selection=(amplitude(mad)>0);/; \
	   s/op\.3\.expression=.*;/op\.3\.expression=(hlnewb=pb);/; \
	   s/op\.3\.selection=.*;/op\.3\.selection=(amplitude(mad)>0);/; \
	   s/op\.4\.expression=.*;/op\.4\.expression=(hlnewc=pc);/; \
	   s/op\.4\.selection=.*;/op\.4\.selection=(amplitude(mad)>0);/; \
	   s/op\.5\.expression=.*;/op\.5\.expression=(hlnewd=pd);/; \
	   s/op\.5\.selection=.*;/op\.5\.selection=(amplitude(mad)>0);/; \
	   s/analysis\.1\.property=\".*\"/analysis\.1\.property=\"average\"/; \
	   s/analysis\.1\.expression=.*;/analysis\.1\.expression=(mad);/; \
	   s/analysis\.1\.selection=.*;/analysis\.1\.selection=(d>4);/; \
	   s/analysis\.2\.property=\".*\"/analysis\.2\.property=\"max\"/; \
	   s/analysis\.2\.expression=.*;/analysis\.2\.expression=(hlnewa);/; \
	   s/analysis\.2\.selection=.*;/analysis\.2\.selection=(d>3);/; \
	   s/analysis\.3\.property=\".*\"/analysis\.3\.property=\"element\"/; \
	   s/analysis\.3\.expression=.*;/analysis\.3\.expression=(hlnewa);/; \
	   s/analysis\.3\.selection=.*;/analysis\.3\.selection=(d>3);/; \
	   s/{===>} reflection_outfile=\".*\"/{===>} reflection_outfile=\"manipulate_1.hkl\"/; \
	   s/{===>} analyse_outfile=\".*\"/{===>} analyse_outfile=\"manipulate_1.list\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

merge.sed:
	@ sed -e \
	  "s/new_name_1_1=\".*\"/new_name_1_1=\"\"/; \
	   s/new_name_2_1=\".*\"/new_name_1_2=\"\"/; \
	   s/{===>} reflection_outfile=\".*\"/{===>} reflection_outfile=\"merge_1.hkl\"/"  \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/reflection_infile_1=\".*\"/reflection_infile_1=\"mbp_refine.hkl\"/; \
	   s/reflection_infile_2=\".*\"/reflection_infile_2=\"\"/; \
	   s/reflection_infile_3=\".*\"/reflection_infile_3=\"\"/; \
	   s/reflection_infile_4=\".*\"/reflection_infile_4=\"\"/; \
	   s/old_name_1_1=\"fobs\"/old_name_1_1=\"pa\"/; \
	   s/old_name_2_1=\"sigma\"/old_name_2_1=\"pb\"/; \
	   s/old_name_3_1=\"\"/old_name_3_1=\"pc\"/; \
	   s/old_name_4_1=\"\"/old_name_4_1=\"pd\"/; \
	   s/old_name_5_1=\"\"/old_name_5_1=\"mad\"/; \
	   s/new_name_1_1=\"f_w1\"/new_name_1_1=\"pa\"/; \
	   s/new_name_2_1=\"s_w1\"/new_name_2_1=\"pb\"/; \
	   s/new_name_5_1=\"\"/new_name_5_1=\"mad1\"/; \
	   s/new_type_1_1=\"real\"/new_type_1_1=\"real\"/; \
	   s/new_type_5_1=\"\"/new_type_5_1=\"complex\"/; \
	   s/{===>} reflection_outfile=\".*\"/{===>} reflection_outfile=\"merge_2.hkl\"/"  \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)

model_phase.sed:

scale.sed:
	@ sed -e \
	  "s/{===>} list_outfile=\".*\"/{===>} list_outfile=\"OUTPUT\"/"  \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

matthews_coef.sed:

transform_map.sed:

average_friedels.sed:

model_fcalc.sed:

flip_friedels.sed:

average_map.sed:

mask_map.sed:

em_map_to_hkl.sed:


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

FIRSTFILES=sdb_manipulate

FILES=density_modify \
      mad_bijvoet_ave \
      optimize_ncsop \
      solvent_mask \
      ir_phase \
      mad_phase \
      cns_to_sdb \
      delete_sites \
      generate_sdb \
      pdb_to_sdb \
      sdb_split \
      sdb_to_pdb \
      sdb_to_sdb \
      shift_sites \
      flip_sites

LINKS=../../xtal_data \
      ../../auxiliary \
      ..

include ../../Makefile.include

density_modify.sed:
	@ sed -e \
	  "s/{===>} sg=\".*\"/{===>} sg=\"P4(1)2(1)2\"/; \
	   s/{===>} a=.*/{===>} a=101.4\;/; \
	   s/{===>} b=.*/{===>} b=101.4\;/; \
	   s/{===>} c=.*/{===>} c=199.5\;/; \
	   s/{===>} alpha=.*/{===>} alpha=90\;/; \
	   s/{===>} beta=.*/{===>} beta=90\;/; \
	   s/{===>} gamma=.*/{===>} gamma=90\;/; \
	   s/{===>} reflection_infile_1=\".*\"/{===>} reflection_infile_1=\"eg1_abcd.cv\"/; \
	   s/{===>} high_res=.*/{===>} high_res=4\;/; \
	   s/{===>} averaging=false/{===>} averaging=true/; \
	   s/{===>} ncs_mask_infile=\".*\"/{===>} ncs_mask_infile=\"eg1.mask\"/; \
	   s/{===>} ncs_infile=\".*\"/{===>} ncs_infile=\"eg1_ncs_strict.def\"/; \
	   s/{===>} solcon=.*/{===>} solcon=0.55\;/; \
	   s/{===>} start_ave_radius=.*/{===>} start_ave_radius=4.5\;/; \
	   s/{===>} finish_ave_radius=.*/{===>} finish_ave_radius=4\;/; \
	   s/{===>} phase_extend=false/{===>} phase_extend=true/; \
	   s/{===>} initial_highres=.*/{===>} initial_highres=4.5\;/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"density_modify_1\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

mad_bijvoet_ave.sed:

optimize_ncsop.sed:

solvent_mask.sed:
	@ sed -e \
	  "s/{===>} mask_method=\".*\"/{===>} mask_method=\"wang-leslie\"/; \
	   s/{===>} ave_radius=4\.0/{===>} ave_radius=8.0/; \
	   s/{===>} mask_outfile=\".*\"/{===>} mask_outfile=\"solvent_mask_wl.omask\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

ir_phase.sed:
	@ sed -e \
	  "s/{===>} sum_method=\".*\"/{===>} sum_method=\"fft\"/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"ir_phase_1\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/{===>} reflection_infile_1=\".*\"/{===>} reflection_infile_1=\"mbp_merge.hkl\"/; \
	   s/{===>} sg=\".*\"/{===>} sg=\"P2(1)2(1)2(1)\"/; \
	   s/{===>} a=.*/{===>} a=65.508\;/; \
	   s/{===>} b=.*/{===>} b=72.216\;/; \
	   s/{===>} c=.*/{===>} c=45.035\;/; \
	   s/{===>} alpha=.*/{===>} alpha=90\;/; \
	   s/{===>} beta=.*/{===>} beta=90\;/; \
	   s/{===>} gamma=.*/{===>} gamma=90\;/; \
	   s/{===>} f_nat=\".*\"/{===>} f_nat=\"f_mod\"/; \
	   s/{===>} s_nat=\".*\"/{===>} s_nat=\"s_w2\"/; \
	   s/{===>} f_1=\".*\"/{===>} f_1=\"f_der\"/; \
	   s/{===>} s_1=\".*\"/{===>} s_1=\"s_w4\"/; \
	   s/{===>} site_segid_1=\".*\"/{===>} site_segid_1=\"MBP\"/; \
	   s/{===>} anom_1=.*/{===>} anom_1=true\;/; \
	   s/{===>} anom_library=\".*\"/{===>} anom_library=\"CNS_XRAYLIB:anom_mo.lib\"/; \
	   s/{===>} high_res_phase=.*/{===>} high_res_phase=3\;/; \
	   s/{===>} high_res_1=.*/{===>} high_res_1=3\;/; \
	   s/{===>} sitedatabase_infile\.1=\".*\"/{===>} sitedatabase_infile\.1=\"mbp_sites.sdb\"/; \
	   s/{===>} sirmir_qstep=.*/{===>} sirmir_qstep=20\;/; \
	   s/{===>} sirmir_fpstep=.*/{===>} sirmir_fpstep=0\;/; \
	   s/{===>} sirmir_fdpstep_ano=.*/{===>} sirmir_fdpstep_ano=20\;/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"ir_phase_2\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)
	@ sed -e \
	  "s/{===>} sum_method=\".*\"/{===>} sum_method=\"fft\"/; \
           s/{===>} reflection_infile_2=\".*\"/{===>} reflection_infile_2=\"ir_phase_1.hkl\"/; \
           s/{===>} priora=\".*\"/{===>} priora=\"pa\"/; \
           s/{===>} priorb=\".*\"/{===>} priorb=\"pb\"/; \
           s/{===>} priorc=\".*\"/{===>} priorc=\"pc\"/; \
           s/{===>} priord=\".*\"/{===>} priord=\"pd\"/; \
           s/{===>} pa_name=\".*\"/{===>} pa_name=\"pa_new\"/; \
           s/{===>} pb_name=\".*\"/{===>} pb_name=\"pb_new\"/; \
           s/{===>} pc_name=\".*\"/{===>} pc_name=\"pc_new\"/; \
           s/{===>} pd_name=\".*\"/{===>} pd_name=\"pd_new\"/; \
           s/{===>} centroid_phase_name=\".*\"/{===>} centroid_phase_name=\"sir_new\"/; \
           s/{===>} fom_name=\".*\"/{===>} fom_name=\"fom_new\"/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"ir_phase_3\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_3.inp)
	@ sed -e \
	  "s/{===>} site_segid_1=\".*\"/{===>} site_segid_1=\"kuo+\"/; \
	   s/{===>} sitedatabase_infile.1=\".*\"/{===>} sitedatabase_infile.1=\"test_sites.sdb\"/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"ir_phase_4\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_4.inp)

mad_phase.sed:
	@ sed -e \
	  "s/{===>} reflection_infile_1=\".*\"/{===>} reflection_infile_1=\"mbp_scale_red.hkl\"/; \
	   s/{===>} high_res=1.7/{===>} high_res=4./; \
	   s/{===>} gradient_map=true/{===>} gradient_map=false/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"mad_phase_1\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/{===>} reflection_infile_1=\".*\"/{===>} reflection_infile_1=\"mad_phase_1.hkl\"/; \
           s/{===>} reflection_infile_2=\".*\"/{===>} reflection_infile_2=\"mbp_scale_red.hkl\"/; \
           s/{===>} priora=\".*\"/{===>} priora=\"pa\"/; \
           s/{===>} priorb=\".*\"/{===>} priorb=\"pb\"/; \
           s/{===>} priorc=\".*\"/{===>} priorc=\"pc\"/; \
           s/{===>} priord=\".*\"/{===>} priord=\"pd\"/; \
           s/{===>} pa_name=\".*\"/{===>} pa_name=\"pa_new\"/; \
           s/{===>} pb_name=\".*\"/{===>} pb_name=\"pb_new\"/; \
           s/{===>} pc_name=\".*\"/{===>} pc_name=\"pc_new\"/; \
           s/{===>} pd_name=\".*\"/{===>} pd_name=\"pd_new\"/; \
           s/{===>} centroid_phase_name=\".*\"/{===>} centroid_phase_name=\"mad_new\"/; \
           s/{===>} fom_name=\".*\"/{===>} fom_name=\"fom_new\"/; \
	   s/{===>} high_res=1.7/{===>} high_res=4./; \
	   s/{===>} gradient_map=true/{===>} gradient_map=false/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"mad_phase_2\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)

sdb_manipulate.sed:

cns_to_sdb.sed:

delete_sites.sed:

generate_sdb.sed:

pdb_to_sdb.sed:

sdb_split.sed:

sdb_to_pdb.sed:

sdb_to_sdb.sed:

shift_sites.sed:

flip_sites.sed:


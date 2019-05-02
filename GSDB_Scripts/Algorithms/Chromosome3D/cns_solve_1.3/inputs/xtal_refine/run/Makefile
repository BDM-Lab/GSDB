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

FIRSTFILES=water_pick

FILES=anneal \
      bdomain \
      bgroup \
      bindividual \
      composite_omit_map \
      fo-fo_map \
      fp_fdp_group \
      minimize \
      map_cover \
      model_map \
      model_stats \
      ncs_average_map \
      optimize_rweight \
      optimize_wa \
      qgroup \
      qindividual \
      rigid \
      sa_omit_map \
      shift_solvent \
      refine \
      water_delete \
      refine_den \
      gradient_map

LINKS=../../xtal_data \
      ../../auxiliary \
      ..

include ../../Makefile.include

anneal.sed:
	@ sed -e \
	  "s/{===>} temperature=[0-9]*/{===>} temperature=300/; \
	   s/{===>} num_trials=[0-9]*/{===>} num_trials=1/; \
	   s/{===>} init_min=[0-9]*/{===>} init_min=10/; \
	   s/{===>} final_min=[0-9]*/{===>} final_min=10/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"anneal_con\"/; \
	   s/{===>} md_scheme=\".*\"/{===>} md_scheme=\"constant\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/{===>} temperature=[0-9]*/{===>} temperature=300/; \
	   s/{===>} num_trials=[0-9]*/{===>} num_trials=1/; \
	   s/{===>} init_min=[0-9]*/{===>} init_min=10/; \
	   s/{===>} final_min=[0-9]*/{===>} final_min=10/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"anneal_cart\"/; \
	   s/{===>} md_type=\".*\"/{===>} md_type=\"cartesian\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)
	@ sed -e \
	  "s/{===>} output_root=\".*\"/{===>} output_root=\"anneal_den\"/; \
	   s/{===>} use_den=false/{===>} use_den=true/" \
	   $(@:.sed=.inp) > $(@:.sed=_3.inp)

bdomain.sed:

bgroup.sed:

bindividual.sed:

composite_omit_map.sed:
	@ sed -e \
	  "s/{===>} output_root=\".*\"/{===>} output_root=\"composite_omit_map_den\"/; \
	   s/{===>} use_den=false/{===>} use_den=true/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)

fo-fo_map.sed:

fp_fdp_group.sed:

minimize.sed:
	@ sed -e \
	  "s/{===>} bulk_sol=true/{===>} bulk_sol=false/; \
	   s/{===>} coordinate_outfile=\".*\"/{===>} coordinate_outfile=\"minimize_nobulk\.pdb\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/{===>} coordinate_outfile=\".*\"/{===>} coordinate_outfile=\"minimize_den.pdb\"/; \
	   s/{===>} use_den=false/{===>} use_den=true/" \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)

map_cover.sed:

model_map.sed:
	@ sed -e \
	  "s/{===>} output_root=\".*\"/{===>} output_root=\"model_map_std\"/; \
	   s/{===>} peak_search=true/{===>} peak_search=false/; \
	   s/{===>} write_coeff=false/{===>} write_coeff=true/; \
	   s/{===>} map_type=\".*\"/{===>} map_type=\"unweighted\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/{===>} output_root=\".*\"/{===>} output_root=\"model_map_2\"/; \
	   s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"amy_ensemble\"/;\
	   s/{===>} coordinate_ensemble_count=1/{===>} coordinate_ensemble_count=3/" \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)

model_stats.sed:

ncs_average_map.sed:
	@ sed -e \
	  "s/{===>} temperature=[0-9]*/{===>} temperature=300/; \
	   s/{===>} cool_rate=[0-9]*/{===>} cool_rate=50/; \
	   s/{===>} obs_pa=\".*\"/{===>} obs_pa=\"pa\"/; \
	   s/{===>} obs_pb=\".*\"/{===>} obs_pb=\"pb\"/; \
	   s/{===>} obs_pc=\".*\"/{===>} obs_pc=\"pc\"/; \
	   s/{===>} obs_pd=\".*\"/{===>} obs_pd=\"pd\"/; \
	   s/{===>} mini_steps=[0-9]*/{===>} mini_steps=10/; \
	   s/{===>} reftarget=\".*\"/{===>} reftarget=\"mlhl\"/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"ncs_comb_map\"/; \
	   s/{===>} map_type=\".*\"/{===>} map_type=\"combined\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/{===>} ncs_infile=\".*\"/{===>} ncs_infile=eg1_ncs_strict.def/; \
	   s/{===>} mask_infile_1=\".*\"/{===>} mask_infile_1=\"\"/; \
	   s/{===>} structure_infile=\".*\"/{===>} structure_infile=\"eg1.mtf\"/; \
	   s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"eg1.pdb\"/; \
	   s/{===>} ncs_group_1=(none)/{===>} ncs_group_1=(all)/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"ncs_strict_map\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)

optimize_rweight.sed:

optimize_wa.sed:

qgroup.sed:

qindividual.sed:

refine.sed:
	@ sed -e \
	  "s/{===>} bfactor_mode=\".*\"/{===>} bfactor_mode=\"group\"/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"refine_1\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_1.inp)
	@ sed -e \
	  "s/{===>} output_root=\".*\"/{===>} output_root=\"refine_2\"/; \
	   s/{===>} write_map=true/{===>} write_map=false/; \
	   s/{===>} full_vdw_cycles=3/{===>} full_vdw_cycles=2/; \
	   s/{===>} den_infile=\".*\"/{===>} den_infile=\"den_make.den\"/; \
	   s/{===>} use_den=false/{===>} use_den=true/" \
	   $(@:.sed=.inp) > $(@:.sed=_2.inp)
	@ sed -e \
	  "s/{===>} bfactor_mode=\".*\"/{===>} bfactor_mode=\"domain\"/; \
	   s/{===>} reset_b=-1/{===>} reset_b=30/; \
	   s/{===>} bd_group_1=(none)/{===>} bd_group_1=(resid 1:20)/; \
	   s/{===>} bd_group_2=(none)/{===>} bd_group_2=(resid 30:40)/; \
	   s/{===>} write_map=true/{===>} write_map=false/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"refine_3\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_3.inp)
	@ sed -e \
	  "s/{===>} minimize_nstep=[0-9]*/{===>} minimize_nstep=0/; \
	   s/{===>} atom_rigid_1=(none)/{===>} atom_rigid_1=(resid 1:34)/; \
	   s/{===>} atom_rigid_2=(none)/{===>} atom_rigid_2=(resid 36:74)/; \
	   s/{===>} write_map=true/{===>} write_map=false/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"refine_4\"/" \
	   $(@:.sed=.inp) > $(@:.sed=_4.inp)
	@ sed -e \
	  "s/{===>} structure_infile=\".*\"/{===>} structure_infile=\"eg1_dimer.mtf\"/; \
	   s/{===>} coordinate_infile=\".*\"/{===>} coordinate_infile=\"eg1_dimer.pdb\"/; \
	   s/{===>} reflection_infile_1=\".*\"/{===>} reflection_infile_1=\"eg1_abcd.cv\"/; \
	   s/{===>} parameter_infile_5=\".*\"/{===>} parameter_infile_5=\"CNS_TOPPAR:carbohydrate.param\"/; \
	   s/{===>} ncs_infile=\".*\"/{===>} ncs_infile=\"ncs.def\"/; \
	   s/{===>} num_cycles=[0-9]*/{===>} num_cycles=2/; \
	   s/{===>} sg=\".*\"/{===>} sg=\"P4(1)2(1)2\"/; \
	   s/{===>} a=[0-9,.]*/{===>} a=101.4/; \
	   s/{===>} b=[0-9,.]*/{===>} b=101.4/; \
	   s/{===>} c=[0-9,.]*/{===>} c=199.5/; \
	   s/{===>} high_res=[0-9,.]*/{===>} high_res=4.0/; \
	   s/{===>} obs_pa=\".*\"/{===>} obs_pa=\"PA\"/; \
	   s/{===>} obs_pb=\".*\"/{===>} obs_pb=\"PB\"/; \
	   s/{===>} obs_pc=\".*\"/{===>} obs_pc=\"PC\"/; \
	   s/{===>} obs_pd=\".*\"/{===>} obs_pd=\"PD\"/; \
	   s/{===>} reftarget=\".*\"/{===>} reftarget=\"mlhl\"/; \
	   s/{===>} write_map=true/{===>} write_map=false/; \
	   s/{===>} anneal=true/{===>} anneal=false/; \
	   s/{===>} output_root=\".*\"/{===>} output_root=\"refine_5\"/ " \
	   $(@:.sed=.inp) > $(@:.sed=_5.inp) 
	@ sed -e \
	  "s/{===>} output_root=\".*\"/{===>} output_root=\"refine_6\"/; \
	   s/{===>} write_map=true/{===>} write_map=false/; \
	   s/{===>} den_infile=\".*\"/{===>} den_infile=\"\"/; \
	   s/{===>} full_vdw_cycles=3/{===>} full_vdw_cycles=2/; \
	   s/{===>} use_den=false/{===>} use_den=true/" \
	   $(@:.sed=.inp) > $(@:.sed=_6.inp)
      

rigid.sed:

sa_omit_map.sed:

shift_solvent.sed:

water_pick.sed:

water_delete.sed:

refine_den.sed:

gradient_map.sed:


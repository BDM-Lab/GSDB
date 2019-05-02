from modeller import *
env = environ()
aln = alignment(env)
mdl = model(env, file='1VGY.pdb', model_segment=('FIRST:A','LAST:A'))
aln.append_model(mdl, align_codes='1VGYA', atom_files='1VGY.pdb')
aln.append(file='target_sequence.pir', align_codes='hp33')
aln.salign(local_alignment=True, rr_file='${LIB}/blosum62.sim.mat',
           gap_penalties_1d=(-600, -600), 
           output='',
           align_block=15,   # no. of seqs. in first MSA
           align_what='PROFILE',
           alignment_type='PAIRWISE',
           comparison_type='PSSM',  # or 'MAT' (Caution: Method NOT benchmarked
                                    # for 'MAT')
           similarity_flag=True,    # The score matrix is not rescaled
           substitution=True,       # The BLOSUM62 substitution values are
                                    # multiplied to the corr. coef.
           #write_weights=True,
           #output_weights_file='test.mtx', # optional, to write weight matrix
           smooth_prof_weight=10.0) # For mixing data with priors

aln.edit(edit_align_codes='hp33', base_align_codes='rest',min_base_entries=1, overhang=0)
aln.write(file='hp33.ali', alignment_format='PIR')
aln.write(file='hp33.pap', alignment_format='PAP')

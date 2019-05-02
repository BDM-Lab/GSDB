from modeller import *

env = environ()
aln = alignment(env)
for (pdb, chain) in (('1cg2', 'A'), ('3ct9', 'A'), ('2f7v', 'A'),
                     ('3gb0', 'A'), ('3isz', 'A'), ('3pfo', 'A'),
                     ('2rb7', 'A'), ('1vgy', 'A') ):
    m = model(env, file=pdb, model_segment=('FIRST:'+chain, 'LAST:'+chain))
    aln.append_model(m, atom_files=pdb, align_codes=pdb+chain)
aln.malign()
aln.malign3d()
aln.compare_structures()
aln.id_table(matrix_file='family.mat')
env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)
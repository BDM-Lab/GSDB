#!/bin/csh
/bin/echo "CNS refinement job with seed $argv[1]"
sed -e 's/{===>} seed=[0-9]*/{===>} seed='$argv[1]'/; \
        s/{===>} den_gamma=[0-9,.]*/{===>} den_gamma='$argv[2]'/; \
        s/{===>} den_scale=[0-9,.]*/{===>} den_scale='$argv[3]'/; \
        s/{===>} minimize_nstep=[0-9,.,-]*/{===>} minimize_nstep=0/; \
        s/{===>} output_root=\".*\"/{===>} output_root=\"refine_den_'$argv[1]'_'$argv[2]'_'$argv[3]'\"/;'  \
        refine_den.inp > refine_den_$argv[1]_$argv[2]_$argv[3].inp

#
# note: change the path for cns_solve_env according to your computing environment
#

source /usr/local/cns/cns_solve_env_local

cns < refine_den_$argv[1]_$argv[2]_$argv[3].inp > refine_den_$argv[1]_$argv[2]_$argv[3].out

rm *.map*
gzip *.out



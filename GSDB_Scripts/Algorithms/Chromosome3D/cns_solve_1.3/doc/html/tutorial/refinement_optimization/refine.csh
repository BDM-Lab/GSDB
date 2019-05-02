#!/bin/csh
/bin/echo "CNS refinement job with seed $argv[1], wa=$argv[2], rweight=$argv[3]"
sed -e 's/{===>} seed=[0-9]*/{===>} seed='$argv[1]'/; \
        s/{===>} wa=[0-9,.,-]*/{===>} wa='$argv[2]'/; \
        s/{===>} rweight=[0-9,.,-]*/{===>} rweight='$argv[3]'/; \
        s/{===>} output_root=\".*\"/{===>} output_root=\"refine_'$argv[1]'_'$argv[2]'_'$argv[3]'\"/;'  \
        refine.inp > refine_$argv[1]_$argv[2]_$argv[3].inp

#
# note: change the path for cns_solve_env according to your computing environment
#
if ( -e /Users/Shared/cns_development/cns_solve_env) source /Users/Shared/cns_development/cns_solve_env

cns < refine_$argv[1]_$argv[2]_$argv[3].inp > refine_$argv[1]_$argv[2]_$argv[3].out 

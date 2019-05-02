#!/bin/csh
/bin/echo "CNS refinement job with seed $argv[1]"
sed -e 's/{===>} seed=[0-9]*/{===>} seed='$argv[1]'/; \
        s/{===>} output_root=\".*\"/{===>} output_root=\"refine_'$argv[1]'\"/;'  \
        refine.inp > refine_$argv[1].inp

#
# note: change the path for cns_solve_env according to your computing environment
#
if ( -e /Users/Shared/cns_development/cns_solve_env) source /Users/Shared/cns_development/cns_solve_env

cns < refine_$argv[1].inp > refine_$argv[1].out 

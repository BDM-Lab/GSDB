#!/bin/bash                                       
echo "starting cns.."                           
touch iam.running                                 
# CNS-CONFIGURATION                               
source /storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/Chromosome3D/cns_solve_1.3/cns_solve_env.sh                
export KMP_AFFINITY=none                          
/storage/htc/bdm/tosin/GSDB_Scripts/Algorithms/Chromosome3D/cns_solve_1.3/intel-x86_64bit-linux/bin/cns_solve < dgsa.inp > /dev/null 
if [ -f "chr18_1mb_matrix_20.pdb" ]; then      
   rm iam.running                                 
   echo "trial structures written."             
   rm *embed*                                     
   exit                                           
fi                                                
if [ -f "chr18_1mb_matrixa_20.pdb" ]; then 
   rm iam.running                                 
   echo "accepted structures written."          
   rm *embed*                                     
   exit                                           
fi                                                
echo "ERROR! Final structures not found!"       
echo "CNS FAILED!"                              
mv iam.running iam.failed                         

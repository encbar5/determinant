#!/bin/bash
# for each number of procs 
 for i in 1 2 4 8 16 32; do
            if [[ $i -lt 20 ]]; then
                ppn=$i
            else
                ppn=20
            fi
            #replace parameters in template sbatch file
            sed -e "s/\#SBATCH --ntasks=.*/\#SBATCH --ntasks=$i/;s/\#SBATCH --ntasks-per-node=.*/\#SBATCH --ntasks-per-node=$ppn/;s/^N_PROC=.*/N_PROC=$i/" $1 > $1.out
            #submit job
            sbatch $1.out
 done
rm $1.out

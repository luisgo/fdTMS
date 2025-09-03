#!/bin/bash -l
#SBATCH --array=0-2%100
#SBATCH --job-name=coilopt
#SBATCH -p wmglab
#SBATCH --cpus-per-task=20
#SBATCH --mem=650000
#SBATCH --time=40:00:00
#SBATCH --mail-user=ljg24@duke.edu
#SBATCH --export=ALL
echo "I ran on:" "${SLURM_ARRAY_TASK_ID}"
#
cd /work/ljg24/loopstar/Arraycode/Simtools/keep
# 
i=$((${SLURM_ARRAY_TASK_ID}))
export OMP_STACKSIZE=4096M     # omp default stack size is too small
ulimit -c 0                    # no core dumps
ulimit -s unlimited            # unlimited stack

/usr/bin/time -v /hpchome/apps/rhel7/matlabR2017b/bin/matlab -nodesktop -nodisplay -r "cd /work/ljg24/BEM_lin2;testcode2($i);exit;"



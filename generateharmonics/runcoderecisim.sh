#!/bin/bash -l
#SBATCH --array=1%100
#SBATCH --job-name=coilopt
#SBATCH -p wmglab
#SBATCH --cpus-per-task=20
#SBATCH --mem=340000
#SBATCH --time=40:00:00
#SBATCH --mail-user=ljg24@duke.edu
#SBATCH --export=ALL
echo "I ran on:" "${SLURM_ARRAY_TASK_ID}"
#
cd /work/ljg24/coiloptcodes/generateharmonics
# 


ulimit -c 0                    # no core dumps
ulimit -s unlimited            # unlimited stack
/hpchome/apps/rhel7/matlabR2017b/bin/matlab -nodesktop -nodisplay -r "cd /work/ljg24/coiloptcodes/generateharmonics;generateharmonicsinputmesh('harmcorrwirethonelayer1126.mat',1);exit;"
#/hpchome/apps/rhel7/matlabR2017b/bin/matlab -nodesktop -nodisplay -r "cd /work/ljg24/coiloptcodes/generateharmonics;generateharmonicsinputmesh('harmcorrwirethtwolayer.mat',2);exit;"

#!/bin/bash
#SBATCH -e /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.err
#SBATCH -o /mnt/lustre/scratch/nlsas/home/csic/eye/lmf/logs/imputation_GLIMPSE/slurm-%j.out

# load the docker module
module load udocker

# run the GLIMPSE image
udocker run --hostenv --hostauth --user=cseyelmf -v /mnt -v /tmp -v=$HOME --workdir=$HOME glimpse:v1.1.1-c27e90d_20210521 /bin/bash ${1}
  
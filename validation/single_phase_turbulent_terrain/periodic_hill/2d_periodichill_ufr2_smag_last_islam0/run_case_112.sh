#!/bin/bash
#SBATCH -A FY260117
#SBATCH -t 48:00:00
#SBATCH --qos=normal
#SBATCH --reservation=flight-cldera

#SBATCH -o sgf_converging_log_%j.out
#SBATCH -J sgf_converging
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ndeveld@sandia.gov
#SBATCH -N 2

source /projects/kynema/builds/default_env_25.sh

# Setup run
nodes=$SLURM_JOB_NUM_NODES
rpn=$RANKS_PER_NODE
ranks=$(( $rpn*$nodes ))

umask 002

srun -N $nodes -n $ranks kynema_sgf terrain_periodic_hill.inp

chown $USER:wg-WindHFM .
chmod g+s .
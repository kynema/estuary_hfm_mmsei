#!/bin/bash
#SBATCH -A FY260117
#SBATCH -t 48:00:00
#SBATCH --qos=normal
#SBATCH --reservation=flight-cldera

#SBATCH -o ugf_ph_log_%j.out
#SBATCH -J ugf_ph
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ndeveld@sandia.gov
#SBATCH -N 3

umask 002
module purge 

module load oneapi/2025.3.0
module load openmpi-oneapi/4.1
module load aue/python/3.12.4
module load gnu/14.2.1

module use /projects/kynema/development/iddesdev_2026-07-10/kynema-manager/spack/share/spack/modules
export RANKS_PER_NODE=112
module load ugf-367f186/kynema-ugf/main-intel-oneapi-compilers-2025.3.0-lqdasxe



# Setup run
nodes=$SLURM_JOB_NUM_NODES
rpn=112
ranks=$(( $rpn * $nodes ))

umask 002

srun -N $nodes -n $ranks kynema_ugf -i periodicHill_smag.yaml

chown $USER:wg-WindHFM .
chmod g+s .

#!/bin/bash
#SBATCH --account=hfm
#SBATCH --nodes=15
#SBATCH --time=08:00:00
#SBATCH -p hbw
#SBATCH --job-name=ros_test_script
#SBATCH --mail-user=mkuhn@nlr.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=%x.%j.out  # %j will be replaced with the job ID

export FI_MR_CACHE_MONITOR=memhooks
export FI_CXI_RX_MATCH_MODE=software
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export MPICH_OFI_NIC_POLICY=NUMA

export KYNEMA_MANAGER=/scratch/mkuhn/kynema-manager
source ${KYNEMA_MANAGER}/start.sh && spack-start
spack env activate -d ${KYNEMA_MANAGER}/environments/env_kynema-sgf
spack install
spack load kynema-sgf
which kynema_sgf

cd /scratch/mkuhn/
cd estuary_flows
cd script_testing
pwd

cp /home/mkuhn/testruns_data/estuary_flows/script_testing/*inp .

srun -N 15 -n 1260 --ntasks-per-node=84 --distribution=block:block --cpu_bind=rank_ldom --kill-on-bad-exit=1 kynema_sgf rosario.inp > rosario.log


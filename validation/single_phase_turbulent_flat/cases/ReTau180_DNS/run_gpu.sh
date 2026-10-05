#!/bin/bash
#SBATCH --job-name=Avg-Re180
#SBATCH -d afterany:14957221
#SBATCH --account=hfm
#SBATCH --exclusive
##SBATCH --partition=debug
##SBATCH --qos=high
#SBATCH --time=08:00:00
#SBATCH --nodes=2
#SBATCH --gpus-per-node=4
#SBATCH --gpu-bind=closest
## Always use one task per GPU (4 gpus per node -> 4 tasks per node)
#SBATCH --ntasks-per-node=4
#SBATCH --output=output.%x_%j.out
#SBATCH --error=error.%x_%j.err
#SBATCH --mail-user=david.montgomery@nlr.gov
#SBATCH --mail-type=ALL

# Load the spack environment
export KYNEMA_MANAGER=/scratch/dmontgo2/kynema-manager
source ${KYNEMA_MANAGER}/start.sh && spack-start
cd ${KYNEMA_MANAGER}/environments/env_kynema_sgf_gpu
quick-activate .
spack load kynema-sgf

# Go to the directory from which our job was launched
cd $SLURM_SUBMIT_DIR

# Max simulation time in hours (should be slightly less than  #SBATCH --time)
max_sim_time=7.9

# Input file, executable, and extra arguments for the job
input_file="turbulent-flat-re-180-averaging.inp"
exec="kynema_sgf"
data_dir=$SLURM_SUBMIT_DIR
latest_chk=$(ls -1v "${data_dir}" | grep "chk" | tail -1)
if [ -n "$latest_chk" ]; then
     extra_args="io.restart_file=${latest_chk}"
else
     echo "No checkpoint files found in ${data_dir}. Starting from initial conditions."
     extra_args=""
fi
extra_args+=" time.max_wall_time=${max_sim_time}"

# Display job info (ntasks is equivalent to number of GPUs here)
echo "Number of Nodes = " $SLURM_JOB_NUM_NODES
echo " Number of GPUS = " $SLURM_NTASKS
echo "Number of Cores = " $SLURM_NTASKS
echo "input file: ${input_file}"
echo "executable: ${exec}"
echo "extra args: ${extra_args}"

# Load modules for GPU job
module purge;
module load PrgEnv-gnu;
module load cuda;
module load craype-x86-milan;

# Run the job
echo "running job"
srun ${exec} ${input_file} ${extra_args}
echo "job has finished" 

#!/bin/bash
set -euo pipefail

# Allocation for job
allocation="hfm"
user_email="david.montgomery@nlr.gov"

# Case details
exe="kynema_sgf"
kynema_manager_dir="/scratch/dmontgo2/kynema-manager"
kynema_env_dir="/scratch/dmontgo2/kynema-manager/environments/env_kynema_sgf"

case_name="max_lev3"
case_dir="/scratch/mkuhn/estuary_flows/milestone/${case_name}"
input_file="rosario_${case_name}.inp"
output_dir="/projects/hfm/dmontgomery/estuary_hfm_mmsei/tools/scalingAnalysis"
tag="${case_name}" # added to job name to distinguish from other runs
max_steps=20010

# Find the latest checkpoint file in the data directory
latest_chk="chk20000"
extra_args="io.restart_file=${case_dir}/${latest_chk}"
extra_args+=" time.max_step=${max_steps} time.plot_interval=-1 time.checkpoint_interval=-1"
extra_args+=" io.plot_file=${output_dir}/${case_name}/plt io.check_file=${output_dir}/${case_name}/chk"

# Turn off all sampling
#extra_args+=" incflo.post_processing=line_sampling line_sampling.output_time_interval=-1"

# Array of node counts to scale over
n_nodes_values=(16)
#n_nodes_values=(32 40 48 56 64)
ntasks_per_node=72

# Temp workspace for generated slurm scripts
tmpdir="$(mktemp -d -t slurm_scale_XXXXXXXX)"
trap 'rm -rf "$tmpdir"' EXIT

for n in "${n_nodes_values[@]}"; do
  job_name="scaleCPU_${tag}_${n}"
  slurm_file="${tmpdir}/${job_name}.slurm"
  ntasks=$((n * ntasks_per_node))

  # Generate a fresh slurm script for this configuration
  cat > "$slurm_file" <<EOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --account=${allocation}
#SBATCH --exclusive
#SBATCH --qos=high
#SBATCH -p hbw
#SBATCH --time=00:10:00
#SBATCH --nodes=${n}
#SBATCH --ntasks-per-node=${ntasks_per_node}
#SBATCH --output=output.%x_%j.out
#SBATCH --mail-user=${user_email}
#SBATCH --mail-type=ALL

set -euo pipefail

export FI_MR_CACHE_MONITOR=memhooks
export FI_CXI_RX_MATCH_MODE=software
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export MPICH_OFI_NIC_POLICY=NUMA

# Load the spack environment
cd ${kynema_manager_dir}
set +u
source start.sh && spack-start
set -u
cd ${kynema_env_dir}
quick-activate .
spack load kynema-sgf

# Go to the directory from which our job was launched
cd \$SLURM_SUBMIT_DIR

echo "        Date: \$(date)"  
echo "      Job ID: \$SLURM_JOB_ID"
echo "       Nodes: \$SLURM_NNODES"
echo "      nTasks: \$SLURM_NTASKS"
echo "  nTasksPerN: \$SLURM_TASKS_PER_NODE"
echo "    location: \$(pwd)"
echo "  input_file: $input_file"
echo "    Hostlist: \$SLURM_NODELIST"

echo "Starting at \$(date)"
srun -N ${n} -n ${ntasks} --ntasks-per-node=${ntasks_per_node} --distribution=block:block --cpu_bind=rank_ldom --kill-on-bad-exit=1 "${exe}" "${input_file}" "${extra_args}"
echo "Finished at \$(date)"
EOF

  # Submit it
  sbatch "$slurm_file"
  echo "Submitted ${job_name} (nodes=${n}, ntasks_per_node=${ntasks_per_node})"
done


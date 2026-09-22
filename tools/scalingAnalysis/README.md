# Scaling Analysis Tools

## `scaling_analysis.py`
Computes and plots HPC scaling metrics (speedup, efficiency, time-to-solution,
AUs) from case definitions in `cases/*.yaml`.

Usage:
```bash
python scaling_analysis.py -i cases/scaling_coarse.yaml
```
Optionally override the final simulation time(s) used for time-to-solution/AU
estimates:
```bash
python scaling_analysis.py -i cases/scaling_coarse.yaml --final-time-hours 12.5 50.0
```
Printed tables go to stdout; plots are saved to `figures/`.

## `mesh_stats.py`
Utility for computing/reporting basic mesh statistics.

## `submit_autoscale_CPU.sh` 
Submit a sweep of SLURM jobs at different node/GPU counts to gather the
run-time data used by `scaling_analysis.py`. Edit the variables near the top
of each script (`input_file`, `tag`, `n_nodes_values`/`n_gpus_values`, etc.)
to match your case, then run:
```bash
./submit_autoscale_CPU.sh
```
Each submits one job per node/GPU count in the sweep array.

Add the `Time spent in Evolve()` output from the `scaleCPU_${tag}_${n}.out` file for each run to a `cases/${tag}.yaml` file for use with `scaling_analysis.py`. The `cases/refine_interface_terrain_adcps_maxlev3.yaml` file can be used as an example.
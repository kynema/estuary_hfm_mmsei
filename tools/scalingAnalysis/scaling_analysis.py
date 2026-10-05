import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import yaml
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

"""
This script calculates and displays HPC scaling results. Ensembles are defined
as YAML files and run via:

    python scaling_analysis.py -i cases/scaling_coarse.yaml

Each YAML file must contain a "name" key plus one or more case entries (see
`cases/*.yaml` for examples). An optional top-level "final_time_hours" key
(a single number or a list of numbers) selects the desired simulation end
time(s) used for the time-to-solution/AUs estimates; if omitted, only
speedup/efficiency are computed. `--final-time-hours` on the command line
overrides whatever is in the YAML file.

Figures are always written to a "figures" directory next to this script,
regardless of the current working directory.
"""

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURES_DIR = SCRIPT_DIR / "figures"


def load_case_file(path):
    """Load a YAML file describing a scaling ensemble"""
    with open(path) as f:
        return yaml.safe_load(f)


class HPCScalingResult:
    @staticmethod
    def nodes_from_gpus(gpus, gpus_per_node=4):
        """Compute the number of nodes needed for a given GPU count

        Assumes `gpus_per_node` GPUs are available per node (default 4).
        """
        gpus = np.asarray(gpus)
        return (gpus + gpus_per_node - 1) // gpus_per_node

    def __init__(
        self,
        num_nodes,
        num_steps,
        run_time,
        mode="CPU",  # CPU or GPU
        num_cores_per_node=128,
        gpus=0,
        gpus_per_node=4,
        qos=1,
        cpu_charge_factor=10,
        gpu_charge_factor=100,
        single_node_single_proc_time=None,
        dt=None,
        final_time_hours=None,
    ):

        if isinstance(num_nodes, np.ndarray):
            self.num_nodes = num_nodes
        else:
            self.num_nodes = np.array(num_nodes)
        self.num_cores_per_node = num_cores_per_node
        self.num_cores = self.num_nodes * num_cores_per_node
        self.num_steps = num_steps
        if isinstance(run_time, np.ndarray):
            self.run_time = run_time
        else:
            self.run_time = np.array(run_time)
        self.mode = mode.upper()
        if self.mode == "CPU":
            self.gpus = 0
            self.gpus_per_node = 0
            self.cpu_charge_factor = cpu_charge_factor
        else:
            self.gpus = np.array(gpus)
            self.gpus_per_node = gpus_per_node
            # Default charge factors are for Kestrel
            # https://www.nrel.gov/hpc/system-resource-allocation-unit
            # The GPU charge assumes non-exclusive GPU usage for single GPU nodes
            # For jobs with more than 1 GPU node, we assume 4 GPUs per node
            self.gpu_charge_fraction = np.clip(self.gpus / self.gpus_per_node, 0, 1)
            self.gpu_charge_factor = gpu_charge_factor

        self.qos = qos
        self.single_node_single_proc_time = single_node_single_proc_time
        self.dt = dt
        self.final_time_hours = final_time_hours
        self._validate()

    def _validate(self):
        if self.mode not in ["CPU", "GPU"]:
            raise ValueError("Mode must be 'CPU' or 'GPU'")
        if len(self.num_nodes) != len(self.run_time) or len(self.num_cores) != len(
            self.run_time
        ):
            raise ValueError("All input lists must be of the same length")

    def AUs(self, wall_time_hours, num_nodes):
        """Compute AUs provided wall time in hours"""
        if self.mode == "CPU":
            charge_factor = self.cpu_charge_factor
        else:
            charge_factor = self.gpu_charge_fraction * self.gpu_charge_factor
        return wall_time_hours * num_nodes * self.qos * charge_factor

    def compute_scaling(self):
        """Compute speedup and efficiency

        If a single-node, single-processor reference time was provided at
        construction time, speedup/efficiency are computed relative to that
        (ideal) baseline. Otherwise they fall back to being relative to the
        smallest case in this run (self.run_time[0]).
        """
        P = self.num_cores if self.mode == "CPU" else self.gpus

        if self.single_node_single_proc_time is not None:
            baseline_time = self.single_node_single_proc_time
        else:
            baseline_time = self.run_time[0]

        speedup = baseline_time / self.run_time
        time_per_step = self.run_time / self.num_steps
        efficiency = speedup / P
        cost_per_step = self.AUs(time_per_step / 3600, self.num_nodes)

        return speedup, efficiency, time_per_step, cost_per_step

    def compute_time_to_solution(self, time_per_step, cost_per_step):
        """Estimate time to solution (days) and total AUs for the full run

        Requires that `dt` (timestep, seconds) and `final_time_hours` (desired
        simulation end time, hours) were provided at construction time.
        """
        if self.dt is None or self.final_time_hours is None:
            return None, None

        total_steps = (self.final_time_hours * 3600) / self.dt
        time_to_solution_days = time_per_step * total_steps / 86400
        total_aus = cost_per_step * total_steps

        return time_to_solution_days, total_aus

    def to_dataframe(self):
        speedup, efficiency, time_per_step, cost = self.compute_scaling()
        time_to_solution_days, total_aus = self.compute_time_to_solution(
            time_per_step, cost
        )

        if self.mode == "CPU":
            data = {
                "Nodes": self.num_nodes,
                "Cores": self.num_cores,
                "Time Per Step": time_per_step,
                "AUs Per Step": cost,
                "Speedup": speedup,
                "Efficiency": efficiency,
            }
        else:
            data = {
                "Nodes": self.num_nodes,
                "GPUs": self.gpus,
                "Time Per Step": time_per_step,
                "AUs Per Step": cost,
                "Speedup": speedup,
                "Efficiency": efficiency,
            }

        if time_to_solution_days is not None:
            data["Total time (days)"] = time_to_solution_days
            data["Total AUs"] = total_aus

        return pd.DataFrame(data)

    def results(self):
        print("-"*75)
        if self.mode == "CPU":
            print("CPU Scaling Results:")
            print("Number of cores per node:", self.num_cores[0] / self.num_nodes[0])

        else:
            print("GPU Scaling Results:")

        if self.single_node_single_proc_time is not None:
            print(
                "Speedup/Efficiency baseline: single-node, single-processor time =",
                self.single_node_single_proc_time,
            )
        else:
            print(
                "Speedup/Efficiency baseline: smallest case in this run (Nodes =",
                f"{self.num_nodes[0]})",
            )

        df = self.to_dataframe()
        two_decimal_cols = {"Speedup", "Total AUs"}
        formatters = {
            col: (lambda x: f"{x:.2f}") if col in two_decimal_cols else (lambda x: f"{x:.3f}")
            for col in df.columns
            if pd.api.types.is_float_dtype(df[col])
        }
        print(df.to_string(index=False, formatters=formatters))
        print("-"*75)


def run_scaling_analysis(cases, final_time_hours=None):
    """Compute, print, and plot scaling results for a dict of cases

    `cases` must contain a "name" key (used as the plot filename prefix)
    plus one or more case entries with the raw scaling inputs. Each case must
    provide a "mode" key ("CPU" or "GPU") and may provide a
    "single_node_single_proc_time" key. CPU cases need "num_nodes" and
    "num_cores_per_node". GPU cases need "gpus" (and optionally
    "gpus_per_node", default 4); the number of nodes is computed via
    `HPCScalingResult.nodes_from_gpus`. Any non-dict entries (e.g. "name" or
    metadata comments like "_comment") are ignored.
    """
    for label, case in cases.items():
        if not isinstance(case, dict):
            continue
        mode = case["mode"].upper()
        single_node_single_proc_time = case.get("single_node_single_proc_time")
        if mode == "CPU":
            case["result"] = HPCScalingResult(
                case["num_nodes"],
                case["num_steps"],
                case["run_times"],
                mode="CPU",
                num_cores_per_node=case["num_cores_per_node"],
                single_node_single_proc_time=single_node_single_proc_time,
                dt=case["dt"],
                final_time_hours=final_time_hours,
            )
        else:
            gpus_per_node = case.get("gpus_per_node", 4)
            num_nodes = HPCScalingResult.nodes_from_gpus(case["gpus"], gpus_per_node)
            case["result"] = HPCScalingResult(
                num_nodes,
                case["num_steps"],
                case["run_times"],
                mode="GPU",
                gpus=case["gpus"],
                gpus_per_node=gpus_per_node,
                single_node_single_proc_time=single_node_single_proc_time,
                dt=case["dt"],
                final_time_hours=final_time_hours,
            )
        case["result"].results()

    plot_scaling(cases)


def plot_scaling(cases):
    """Create scaling plots vs. Nodes for a dict of {label: case} entries"""

    case_description = cases["name"]
    final_time_hours = cases[list(cases.keys())[1]]['result'].final_time_hours
    dt = cases[list(cases.keys())[1]]['result'].dt

    fig1, (ax_speedup, ax_eff) = plt.subplots(1, 2, figsize=(12, 5))
    fig2, (ax_time, ax_aus) = plt.subplots(1, 2, figsize=(12, 5))

    for label, case in cases.items():
        if not isinstance(case, dict):
            continue
        df = case["result"].to_dataframe()
        marker = "s" if case["mode"].upper() == "GPU" else "o"

        ax_speedup.plot(df["Nodes"], df["Speedup"], marker=marker, label=label)
        ax_eff.plot(df["Nodes"], df["Efficiency"], marker=marker, label=label)

        if "Total time (days)" in df.columns:
            ax_time.plot(df["Nodes"], df["Total time (days)"], marker=marker, label=label)
            ax_aus.plot(df["Nodes"], df["Total AUs"], marker=marker, label=label)

    ax_speedup.set_xlabel("Nodes")
    ax_speedup.set_ylabel("Speedup")
    ax_speedup.set_title("Speedup vs. Nodes")
    ax_speedup.legend()
    ax_speedup.grid(True)

    ax_eff.set_xlabel("Nodes")
    ax_eff.set_ylabel("Efficiency")
    ax_eff.set_title("Efficiency vs. Nodes")
    ax_eff.legend()
    ax_eff.grid(True)

    ax_time.set_xlabel("Nodes")
    ax_time.set_ylabel("Est. Simulation Time (days)")
    ax_time.set_title(f"Est. Simulation Time vs. Nodes ($t_f$ = {final_time_hours} hrs)")
    ax_time.legend()
    ax_time.grid(True)

    ax_aus.set_xlabel("Nodes")
    ax_aus.set_ylabel("Est. AUs to Solution")
    ax_aus.set_title(f"Est. AUs to Solution vs. Nodes ($t_f$ = {final_time_hours} hrs)")
    ax_aus.legend()
    ax_aus.grid(True)

    for ax in (ax_speedup, ax_eff, ax_time, ax_aus):
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    fig1.tight_layout()
    fig2.tight_layout()
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig1.savefig(FIGURES_DIR / f"{case_description}_performance_tf_{final_time_hours}.png")
    fig2.savefig(FIGURES_DIR / f"{case_description}_estimates_tf_{final_time_hours}.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run HPC scaling analysis for an ensemble of cases defined in a YAML file"
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help=(
            "Path to a YAML file describing the scaling cases (see cases/*.yaml). "
            "Relative paths are resolved relative to this script's directory."
        ),
    )
    parser.add_argument(
        "--final-time-hours", type=float, nargs="+", default=None,
        help=(
            "Desired simulation end time(s) in hours, used for time-to-solution/AUs "
            "estimates. Pass multiple values to run the analysis once per value. "
            "Overrides any 'final_time_hours' entry in the YAML file."
        ),
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_absolute():
        input_path = SCRIPT_DIR / input_path

    cases = load_case_file(input_path)

    final_time_hours_values = args.final_time_hours
    if final_time_hours_values is None:
        final_time_hours_values = cases.pop("final_time_hours", None)
        if final_time_hours_values is None:
            final_time_hours_values = [None]
        elif not isinstance(final_time_hours_values, list):
            final_time_hours_values = [final_time_hours_values]
    else:
        cases.pop("final_time_hours", None)

    for final_time_hours in final_time_hours_values:
        print(f"\n\n\nScaling analysis for final_time_hours = {final_time_hours} hrs")
        run_scaling_analysis(cases, final_time_hours=final_time_hours)
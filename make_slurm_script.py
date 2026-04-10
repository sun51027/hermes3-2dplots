#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
import stat
import subprocess

scratch_path="/users/jpm590/scratch"

TEMPLATE = """#!/bin/bash
#
#SBATCH --job-name={job_name}
#SBATCH --nodes={nodes}
#SBATCH --tasks-per-node={cores_per_node}
#SBATCH --time={time}
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/users/jpm590/scratch/slurmlog/%x-%j.log
#SBATCH --error=/users/jpm590/scratch/slurmlog/%x-%j.err
#SBATCH --account=pet-pic-2022
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jpm590@york.ac.uk

set -e
module purge
source /users/jpm590/neutralimit_gcc12_9497476/bout.env

hermes_path=/users/jpm590/neutralimit_gcc12_9497476/hermes-3/build-mc-master/hermes-3
input_base=/users/jpm590/scratch/{input_base}

DRY_RUN="${{DRY_RUN:-${{dryrun:-false}}}}"

if [ "$DRY_RUN" = true ]; then
    if [ ! -f "$input_base/BOUT.inp" ]; then
        echo "[DRY RUN] ERROR: $input_base/BOUT.inp not found." >&2
        exit 1
    fi
    echo "[DRY RUN] BOUT.inp found."
    echo "[DRY RUN] Would run: {runcommand}"
    exit 0
fi

if [ ! -f "$input_base/BOUT.inp" ]; then
    echo "ERROR: $input_base/BOUT.inp not found." >&2
    exit 1
fi

{runcommand}
"""

def main():
    parser = argparse.ArgumentParser(
        description="Generate SLURM script for hermes-3 with restart/dry-run support"
    )
    parser.add_argument("-j",   "--job-name", required=False, help="job's name: YY-MM-DD-2D-MASTU-test", )
    parser.add_argument("-n",   "--nodes", type=int, required=True, help='nodes')
    parser.add_argument("-c",   "--cores_per_node", type=int, required=True, help="number of core per node")
    parser.add_argument("-t",   "--time", required=True, help="running duration")
    parser.add_argument("-i",   "--input-base", required=True, help="input files BOUT.inp folder")
    parser.add_argument("-d",   "--date", required=True, help="YY:MM:DD")
    parser.add_argument("--mkdir", action="store_true", help="Create a new folder and copy restart file from 260217-cdn-46895-nowallpump_limiterOff_1e21 and input file at scratch")
    parser.add_argument("--restart", action="store_true")
    parser.add_argument("--append", action="store_true")
    # parser.add_argument("-task","--ntasks", type=int, required=False, help="(don't use it) number of core for the simulation")
    # parser.add_argument("-c",   "--cpus-per-task", type=int, required=True, help="number of CPU for each task")

    args = parser.parse_args()

    # Enforce time must be explicitly given
#    if not args.time:
#        sys.exit("You need to specify time (e.g. --time 04:00:00)")

    # Default job-name from final folder of input-base if not specified
    if not args.job_name:
        args.job_name = Path(args.input_base).name
        if not args.job_name:
            sys.exit("Could not determine job-name from input-base. Please specify --job-name explicitly.")


    restart_suffix = " restart" if args.restart else ""
    append_suffix = " append" if args.append else ""
    #restart_suffix = " restart append" if args.restart else ""
    total_cores = args.nodes * args.cores_per_node
    script_text = TEMPLATE.format(
        job_name=args.job_name,
        nodes=args.nodes,
        cores_per_node=args.cores_per_node,
        time=args.time,
        input_base=args.input_base,
        restart_suffix=restart_suffix,
        append_suffix=append_suffix,
        runcommand = f'mpirun -np {total_cores} "$hermes_path" -d "$input_base" {restart_suffix} {append_suffix}'
        #ntasks=args.ntasks,
        #cpus_per_task=args.cpus_per_task,
    )

    if args.mkdir:
    # 1. Fixed typo (commends -> commands)
    # 2. Used triple quotes for multi-line formatting
    # 3. Added backslashes (\) to prevent the shell from seeing 
    #    actual newlines/tabs between commands.
        commands = f"""mkdir -p {scratch_path}/{args.input_base} \
            && cp {scratch_path}/260217-cdn-46895-nowallpump_limiterOff_1e21/*restart* {scratch_path}/{args.input_base}/ \
            && cp {scratch_path}/BOUT.inp {scratch_path}/{args.input_base}/"""
        print(commands)
        subprocess.call(commands, shell=True)

    out_dir = Path("/users/jpm590/neutralimit_gcc12_9497476/hermes-3/build-mc-master")
    # out_dir = Path("/users/jpm590/2dspace/hermes-3/build-mc-master")
    out_dir.mkdir(parents=True, exist_ok=True)
    script_path = out_dir / f"job_{args.job_name}.sh"
    script_path.write_text(script_text)
    script_path.chmod(script_path.stat().st_mode | stat.S_IXUSR)


    print(f"✅ Wrote SLURM script: {script_path}")
    print("Examples:")
    print(f"dryrun=true  bash {script_path}")
    print(f"dryrun=false sbatch {script_path}")

if __name__ == "__main__":
    main()


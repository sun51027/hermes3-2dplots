#!/usr/bin/env python3
"""
Apply neutral density (Nd_src) and heat (Pd_src) sources to a Hermes-3 grid file,
then report total integrated quantities before and after.

Usage example:
  python apply_puff_sources.py \
    --old-grid /path/to/old.nc \
    --new-grid /path/to/new.nc \
    --Nd 7.5e21 \
    --Pd 1.5e21
"""

from pathlib import Path
import argparse
import numpy as np
import sys
sys.path.append(r"/users/jpm590/2dspace/post-processing/sdtools")

from hermes3.case_db import *
from hermes3.load import *
from hermes3.named_selections import *
from hermes3.plotting import *
from hermes3.grid_fields import *
from hermes3.accessors import *
from hermes3.utils import *
# ---------- helper functions ----------

def compute_cell_volumes(mesh: Mesh) -> np.ndarray:
    """Toroidal cell volume: dv = dy * dx * J * 2π."""
    return mesh.mesh["dy"] * mesh.mesh["dx"] * mesh.mesh["J"] * (2.0 * np.pi)


def inner_midplane_region(mesh: Mesh):
    """Return two-cell puff region on the inner midplane."""
    j1, j2 = int(mesh.j1_1g), int(mesh.j2_1g)
    imp_a = int((j2 - j1) / 2) + j1 + 1
    imp_b = int((j2 - j1) / 2) + j1
    return (-mesh.MXG - 1, np.r_[imp_a, imp_b])

def outer_midplane_region(mesh: Mesh):
    """Return two-cell puff region on the inner midplane."""
    '''See     sdtools/hermes3/load.py '''

    omp_a = int((mesh["j2_2g"] - mesh["j1_2g"]) / 2) + mesh["j1_2g"]
    omp_b = int((mesh["j2_2g"] - mesh["j1_2g"]) / 2) + mesh["j1_2g"] + 1
    return (-mesh.MXG - 1, np.r_[omp_a, omp_b])

def write_pump_mask(mesh: Mesh) -> None:
    """Mark pump regions (targets and SOL/PFR edges)."""
    is_pump = Field("is_pump", mesh)
    for region in [
        "inner_lower_target", "inner_upper_target",
        "outer_lower_target", "outer_upper_target",
        "sol_edge", "pfr_edge",
    ]:
        is_pump.data[mesh.slices(region)] = 1
    mesh.write_field(is_pump, dtype="Field2D")


def summarize_sources(grid_path: Path, label: str):
    """Return total Nd and Pd sources (integrated over volume)."""
    mesh = Mesh(str(grid_path))
    try:
        dv = compute_cell_volumes(mesh)
        Pd = mesh.mesh["Pd_src"].squeeze()
        Nd = mesh.mesh["Nd_src"].squeeze()
        total_P = float((Pd * dv).sum())
        total_N = float((Nd * dv).sum())
        print("-" * 48)
        print(f"{label}")
        print(f"  Total Pd_src: {total_P:.3e} [W]")
        print(f"  Total Nd_src: {total_N:.3e} [s^-1]")
        return total_P, total_N
    finally:
        try:
            mesh.close()
        except Exception:
            pass
        del mesh


def apply_sources(new_grid_path: Path, total_P: float, total_N: float, puff_loc: str):
    """Apply distributed sources on new grid based on total values."""
    mesh = Mesh(str(new_grid_path))
    try:
        if puff_loc == "omp":
            puff_region = outer_midplane_region(mesh)
        elif puff_loc == "imp":
            puff_region = inner_midplane_region(mesh)

        dv = compute_cell_volumes(mesh)
        region_vol = dv[puff_region]

        # distribute proportionally to volume
        Nd_src = Field("Nd_src", mesh)
        Pd_src = Field("Pd_src", mesh)

        Nd_src.data[puff_region] = total_N * (region_vol / region_vol.sum()) / region_vol
        Pd_src.data[puff_region] = total_P * (region_vol / region_vol.sum()) / region_vol

        mesh.write_field(Nd_src)
        mesh.write_field(Pd_src)
        write_pump_mask(mesh)

        print(f"Applied total Nd_src: {total_N:.3e} [s^-1]")
        print(f"Applied total Pd_src: {total_P:.3e} [W]")
    finally:
        try:
            mesh.close()
        except Exception:
            pass
        del mesh


# ---------- main ----------

def parse_args():
    p = argparse.ArgumentParser(
        description="Impose sources from old grid to new grid for Hermes-3 runs."
    )
    p.add_argument("--old-grid", required=True, type=Path, help="Path to old grid .nc file")
    p.add_argument("--new-grid", required=True, type=Path, help="Path to new grid .nc file")
    p.add_argument("--Nd", type=float, default=None,
                   help="Total neutral influx [s^-1]. If not given, use value from old grid.")
    p.add_argument("--Pd", type=float, default=None,
                   help="Total heat power [W]. If not given, use value from old grid.")
    p.add_argument("--puff", choices=["omp", "imp"],  required=True, type=str,
                   help="puff location inner midplnae (imp) or outer midplane (omp)")
    p.add_argument("--mode", choices=["show", "edit"], default="show",
                   help="'show' to only print summary (default), 'edit' to apply sources.")

    # Optional field coefficients
    p.add_argument("--D-core", type=float, default=4.0)
    p.add_argument("--chi-core", type=float, default=10.0)
    p.add_argument("--D-sol", type=float, default=4.0)
    p.add_argument("--chi-sol", type=float, default=4.0)
    p.add_argument("--Ni-src-core", type=float, default=0.0)
    p.add_argument("--Pi-src-core", type=float, default=0.0)
    p.add_argument("--Pe-src-core", type=float, default=0.0)
    return p.parse_args()


def main():
    args = parse_args()
    # Step 1: always show current totals
    total_P_old, total_N_old = summarize_sources(args.old_grid, "Old grid")

    if args.mode == "show":
        print("\nMode: show (no file modification)")
        print("→ Only showing source summaries.")
        return

    # mode == "edit"
    print("\nMode: edit (fields will be updated)\n")

    # Step 1. Impose fields from old grid to new
    impose_fields(
        str(args.old_grid),
        str(args.new_grid),
        Ni_src_core=args.Ni_src_core,
        Pi_src_core=args.Pi_src_core,
        Pe_src_core=args.Pe_src_core,
        D_core=args.D_core,
        chi_core=args.chi_core,
        D_sol=args.D_sol,
        chi_sol=args.chi_sol,
    )

    # Step 2. Get total Nd, Pd from old grid if not provided
    total_P_old, total_N_old = summarize_sources(args.old_grid, "Old grid")
    total_P = args.Pd if args.Pd is not None else total_P_old
    total_N = args.Nd if args.Nd is not None else total_N_old

    # Step 3. Apply sources on the new grid
    apply_sources(args.new_grid, total_P=total_P, total_N=total_N, puff_loc=args.puff)

    # Step 4. Verify results
    summarize_sources(args.new_grid, "New grid")


if __name__ == "__main__":
    main()


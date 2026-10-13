"""
Check how SOLPS and Hermes-3 treat the target guard cells.

Script version of the "Compare important parameters between SOLPS and Hermes-3"
section in plot_solps_balance.ipynb. For each target it prints:

    A. Total ion flux to the target (SOLPS vs Hermes-3)
    B. Sheath face area, ring by ring
    C. SOLPS guard cell: thickness and Mach number
    D. Hermes-3 face (target) vs SOLPS guard (target): Mach number and values
    E. (commented out) Hermes-3 sheath-extrapolated face value vs SOLPS
    F. Verdict: is the target value sonic in each code, and do they match

Usage (from hermes3-2dplots/):
    python check_guard_cell.py
    python check_guard_cell.py -i 261001_highrsl_fixplsbg_fromsolps_1e21_roomT_gd \
        -s ../SOLPS_origin_1e21 -t outer_lower inner_lower

Naming convention used below:
    "last domain cell"   = the last real (physical) cell before the target
    "second last cell"   = one cell further away from the target
    "guard cell"         = the first cell behind the target (outside the domain)
    "face"               = the boundary between last domain cell and guard cell (= target)
"""

import argparse
import sys

import numpy as np
import pandas as pd

sys.path.append(r"/Users/zero/workspace/phd_work/hms_output/sdtools/")
sys.path.append(r"/Users/zero/workspace/phd_work/hms_output/hermes3-2dplots/pltool")
from code_comparison.solps_pp import *
from hermes3.case_db import *
from hermes3.load import *
from hermes3.named_selections import *
from hermes3.accessors import *
from hermes3.accessors import _select_region   # returns the (x, theta) index slices of a named region
from hermes3.utils import *

from common import adapt_solps_conventions

ELECTRON_CHARGE = 1.602e-19   # [C]
PROTON_MASS     = 1.6726e-27  # [kg]
TARGET_NAMES    = ["inner_lower", "inner_upper", "outer_upper", "outer_lower"]


# ---------------------------------------------------------------------------
# Command line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Check SOLPS vs Hermes-3 target guard cells")
parser.add_argument("-i", "--input", default="261001_highrsl_fixplsbg_fromsolps_1e21_roomT_gd",
                    help="Hermes-3 case ID in hms_output/")
parser.add_argument("-s", "--solps", default="../SOLPS_origin_1e21",
                    help="Path to SOLPS run folder (with balance.nc)")
parser.add_argument("-t", "--targets", nargs="+", default=["outer_lower"], choices=TARGET_NAMES,
                    help="Targets to check")
parser.add_argument("--tol", type=float, default=0.1,
                    help="Tolerance: |M-1| < tol counts as sonic, |ratio-1| < tol counts as matched")
parser.add_argument("--all-rings", action="store_true",
                    help="Print the full ring-by-ring tables (default: first and last 5 rings)")
args = parser.parse_args()
tolerance = args.tol


# ---------------------------------------------------------------------------
# Read SOLPS and Hermes-3 (same as notebook)
# ---------------------------------------------------------------------------
# bal = SOLPS case (same name as in the notebook)
bal = SOLPScase(args.solps)
bal = adapt_solps_conventions(bal)

db = CaseDB(
    case_dir = r"/Users/zero/workspace/phd_work/hms_output/",
    grid_dir = r"/Users/zero/workspace/phd_work/hms_output/gridfile"
)
case = db.load_case_2D(args.input, use_squash = True, verbose = True)
case.extract_2d_tokamak_geometry()

# ds = Hermes-3 dataset at the last time step (same name as in the notebook)
ds = case.ds.isel(t=-1)

hermes_metadata    = ds.metadata
ion_atomic_mass    = ds.options["d+"]["AA"]                         # AA, e.g. 2 for deuterium
recycle_multiplier = ds.options["d+"]["target_recycle_multiplier"]  # e.g. 0.99

# SOLPS poloidal index of the last domain cell and the guard cell, per target.
# NOTE: solps_pp sets psel["inner_lower_target_guard"] = 1 (same as the domain cell);
#       the real guard cell of the inner lower target is index 0.
solps_cell_index = {
    "inner_lower": {"last_domain_cell": bal.psel["inner_lower_target"],
                    "guard_cell":       0},
    "inner_upper": {"last_domain_cell": bal.psel["inner_upper_target"],
                    "guard_cell":       bal.psel["inner_upper_target_guard"]},
    "outer_upper": {"last_domain_cell": bal.psel["outer_upper_target"],
                    "guard_cell":       bal.psel["outer_upper_target_guard"]},
    "outer_lower": {"last_domain_cell": bal.psel["outer_lower_target"],
                    "guard_cell":       bal.psel["outer_lower_target_guard"]},
}

# SOLPS poloidal index of the cell whose west face is the target face
# (same choice as in adapt_solps_conventions)
solps_target_face_index = {
    "inner_lower": bal.psel["inner_lower_target"],
    "inner_upper": bal.psel["inner_upper_target_guard"],
    "outer_upper": bal.psel["outer_upper_target"],
    "outer_lower": bal.psel["outer_lower_target_guard"],
}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
def print_section_header(title, line_character="="):
    """Print a title between two horizontal lines."""
    print("\n" + line_character * 78)
    print(title)
    print(line_character * 78)


def solps_value(variable_name, poloidal_index):
    """
    SOLPS value of one variable along one poloidal index, for all radial rings.

    bal.bal[variable_name] has shape (poloidal, radial).
    The radial range 1:-1 removes the two SOLPS radial guard cells
    (first and last ring), so the result matches the Hermes-3 domain rings.
    """
    return bal.bal[variable_name][poloidal_index, 1:-1]


def compute_mach_number(electron_temperature, ion_temperature, parallel_velocity):
    """Mach number = |V| / cs, with cs = sqrt(e (Te + Ti) / m_i). Temperatures in eV."""
    sound_speed = np.sqrt(ELECTRON_CHARGE * (electron_temperature + ion_temperature)
                          / (ion_atomic_mass * PROTON_MASS))
    return np.abs(parallel_velocity) / sound_speed


def hermes_face_value(last_domain_cell, guard_cell, variable_name):
    """Value Hermes-3 sees at the target face = average of last domain cell and guard cell."""
    return 0.5 * (last_domain_cell[variable_name].values + guard_cell[variable_name].values)


def hermes_parallel_cell_length(cell):
    """Parallel length of a Hermes-3 cell along the field line: dy * sqrt(g_22)."""
    return (cell["dy"] * np.sqrt(cell["g_22"])).values


def compute_sheath_face_area(last_domain_cell, guard_cell):
    """Hermes-3 sheath face area between the last domain cell and the first guard cell."""
    jacobian_over_length = ((last_domain_cell["J"] + guard_cell["J"])
                            / (np.sqrt(last_domain_cell["g_22"]) + np.sqrt(guard_cell["g_22"])))
    average_dx = 0.5 * (last_domain_cell["dx"] + guard_cell["dx"])
    average_dz = 0.5 * (last_domain_cell["dz"] + guard_cell["dz"])
    return (jacobian_over_length * average_dx * average_dz).values


def sheath_extrapolated_face_value(last_domain_cell, second_last_cell, variable_name):
    """Face value after sheath_boundary_simple (mode 1: exponential extrapolation)."""
    value_last        = last_domain_cell[variable_name].values
    value_second_last = second_last_cell[variable_name].values
    return 0.5 * (value_last + value_last**2 / value_second_last)


def count_sonic_rings(mach_number):
    """Return (number of rings with |M - 1| < tolerance, total number of rings)."""
    number_of_sonic_rings = int(np.sum(np.abs(mach_number - 1) < tolerance))
    total_number_of_rings = len(mach_number)
    return number_of_sonic_rings, total_number_of_rings


def print_ratio_summary(label, ratio):
    """One line: mean / min / max of a ratio and how many rings are within tolerance."""
    number_of_matched_rings = int(np.sum(np.abs(ratio - 1) < tolerance))
    if number_of_matched_rings == len(ratio):
        status_symbol = "✅"
    else:
        status_symbol = "❌"
    print(f"  {status_symbol} {label:28s} mean {np.mean(ratio):7.3f}   min {np.min(ratio):7.3f}   "
          f"max {np.max(ratio):7.3f}   within ±{tolerance:.0%}: "
          f"{number_of_matched_rings}/{len(ratio)} rings")


def print_ring_table(table):
    """Print a ring-by-ring table. Long tables show only the first and last 5 rings."""
    with pd.option_context("display.float_format", "{:.3f}".format,
                           "display.width", 200, "display.max_columns", 20):
        if args.all_rings or len(table) <= 10:
            print(table.to_string())
        else:
            first_and_last_rings = pd.concat([table.head(5), table.tail(5)])
            print(first_and_last_rings.to_string())
            print(f"  ... ({len(table)} rings, first and last 5 shown; "
                  f"use --all-rings to print all)")


# ---------------------------------------------------------------------------
# A. Total ion flux to targets
# ---------------------------------------------------------------------------
print_section_header("A. Total ion particle flow to the targets [s^-1]")

number_of_x_guards = hermes_metadata["MXG"]
recycled_flow_per_cell = ds["Sd_target_recycle"] * ds["dv"] / recycle_multiplier
total_flow_all_cells   = float(recycled_flow_per_cell.sum())
total_flow_domain_only = float(recycled_flow_per_cell.isel(
    x=slice(number_of_x_guards, -number_of_x_guards)).sum())

flux_table_rows = []
for target_name in TARGET_NAMES:
    target_cells = ds.hermesm.select_region(f"{target_name}_target").squeeze()
    flow_hermes = float((target_cells["Sd_target_recycle"] * target_cells["dv"]).sum()) / recycle_multiplier
    flow_solps  = bal.target_flux[target_name]
    flux_table_rows.append({"target": target_name, "Hermes3": flow_hermes, "SOLPS": flow_solps,
                            "ratio": flow_hermes / flow_solps})
flux_table_rows.append({"target": "total (domain)", "Hermes3": total_flow_domain_only,
                        "SOLPS": bal.target_flux["total"],
                        "ratio": total_flow_domain_only / bal.target_flux["total"]})

flux_table = pd.DataFrame(flux_table_rows).set_index("target")
print(flux_table.to_string(formatters={"Hermes3": "{:.4e}".format,
                                       "SOLPS":   "{:.4e}".format,
                                       "ratio":   "{:.4f}".format}))
print(f"\n  Flow in x-guard rings (should be ~0): "
      f"{total_flow_all_cells - total_flow_domain_only:.4e}")


# ---------------------------------------------------------------------------
# Per-target checks
# ---------------------------------------------------------------------------
verdicts_per_target = {}

for target_name in args.targets:
    print_section_header(f"TARGET: {target_name}", "#")

    # --- Hermes-3: find the theta index of the last domain cell and the guard cell ---
    # _select_region returns (x_slice, theta_slice); .start gives the first index of the slice
    last_cell_region  = _select_region(ds, f"{target_name}_target")
    guard_cell_region = _select_region(ds, f"{target_name}_target_inner_guards")
    radial_slice_without_guards = last_cell_region[0]

    theta_index_last_cell  = last_cell_region[1].start
    theta_index_guard_cell = guard_cell_region[1].start
    # direction from last cell to guard cell is +1 or -1, so go one step the other way
    step_towards_target = theta_index_guard_cell - theta_index_last_cell
    theta_index_second_last_cell = theta_index_last_cell - step_towards_target

    hermes_last_cell        = ds.isel(theta=theta_index_last_cell,        x=radial_slice_without_guards)
    hermes_guard_cell       = ds.isel(theta=theta_index_guard_cell,       x=radial_slice_without_guards)
    hermes_second_last_cell = ds.isel(theta=theta_index_second_last_cell, x=radial_slice_without_guards)

    # --- SOLPS: poloidal index of the last domain cell and the guard cell ---
    solps_last_index  = solps_cell_index[target_name]["last_domain_cell"]
    solps_guard_index = solps_cell_index[target_name]["guard_cell"]

    # --- Make sure both codes have the same number of radial rings ---
    number_of_rings_hermes = hermes_last_cell.sizes["x"]
    number_of_rings_solps  = solps_value("Te", solps_guard_index).shape[0]
    if number_of_rings_hermes != number_of_rings_solps:
        print(f"❌ Ring count differs: Hermes-3 {number_of_rings_hermes}, "
              f"SOLPS {number_of_rings_solps}. Skipping ring-by-ring comparison for {target_name}.")
        continue
    ring_numbers = pd.Index(np.arange(number_of_rings_hermes), name="ring")

    print(f"Rings: {number_of_rings_hermes}")
    print(f"Hermes-3 theta index: second last = {theta_index_second_last_cell}, "
          f"last = {theta_index_last_cell}, guard = {theta_index_guard_cell}")
    print(f"SOLPS poloidal index: last = {solps_last_index}, guard = {solps_guard_index}")

    # -----------------------------------------------------------------------
    # B. Sheath face area
    # -----------------------------------------------------------------------
    print_section_header("B. Sheath face area (Hermes-3 / SOLPS parallel-projected)", "-")
    solps_face_index = solps_target_face_index[target_name]
    face_area_hermes = compute_sheath_face_area(hermes_last_cell, hermes_guard_cell)
    # gs[..., 0] = poloidal face area; multiply by Bpol/Btot to project it along the field line
    face_area_solps = (bal.bal["gs"][solps_face_index, 1:-1, 0]
                       * np.abs(bal.g["Bpol"][solps_face_index, 1:-1]
                                / bal.g["Btot"][solps_face_index, 1:-1]))
    print_ratio_summary("A_face H3 / SOLPS", face_area_hermes / face_area_solps)

    # -----------------------------------------------------------------------
    # C. SOLPS guard cell
    # -----------------------------------------------------------------------
    print_section_header("C. SOLPS guard cell", "-")
    # hx = poloidal width of a SOLPS cell
    solps_width_ratio_guard_to_last = (bal.g["hx"][solps_guard_index, 1:-1]
                                       / bal.g["hx"][solps_last_index, 1:-1])
    mach_solps_guard = solps_value("M", solps_guard_index)
    mach_solps_last  = solps_value("M", solps_last_index)

    if np.mean(solps_width_ratio_guard_to_last) < 0.1:
        thickness_comment = "thin guard (sits on the target)"
    else:
        thickness_comment = "guard is NOT thin"
    print(f"  Guard / last cell width (hx): mean {np.mean(solps_width_ratio_guard_to_last):.4f}"
          f"  → {thickness_comment}")

    solps_table = pd.DataFrame({"hx guard/last": solps_width_ratio_guard_to_last,
                                "M last":        mach_solps_last,
                                "M guard":       mach_solps_guard},
                               index=ring_numbers)
    print_ring_table(solps_table)

    # -----------------------------------------------------------------------
    # D. Hermes-3 face (target) vs SOLPS guard (target)
    # -----------------------------------------------------------------------
    print_section_header("D. Hermes-3 face (target) vs SOLPS guard (target)", "-")
    hermes_length_ratio_guard_to_last = (hermes_parallel_cell_length(hermes_guard_cell)
                                         / hermes_parallel_cell_length(hermes_last_cell))
    if np.mean(hermes_length_ratio_guard_to_last) < 0.1:
        thickness_comment = "thin guard"
    else:
        thickness_comment = "thick guard: face = (last + guard)/2"
    print(f"  Guard / last cell length (dy*sqrt(g_22)): "
          f"mean {np.mean(hermes_length_ratio_guard_to_last):.3f}  → {thickness_comment}")

    mach_hermes_last = compute_mach_number(hermes_last_cell["Te"].values,
                                           hermes_last_cell["Td+"].values,
                                           hermes_last_cell["Vd+"].values)
    mach_hermes_face = compute_mach_number(
        hermes_face_value(hermes_last_cell, hermes_guard_cell, "Te"),
        hermes_face_value(hermes_last_cell, hermes_guard_cell, "Td+"),
        hermes_face_value(hermes_last_cell, hermes_guard_cell, "Vd+"))

    hermes_table = pd.DataFrame({"dl guard/last":   hermes_length_ratio_guard_to_last,
                                 "M last (H3)":     mach_hermes_last,
                                 "M face (H3)":     mach_hermes_face,
                                 "M guard (SOLPS)": mach_solps_guard,
                                 "M face/guard":    mach_hermes_face / mach_solps_guard},
                                index=ring_numbers)
    print_ring_table(hermes_table)

    print("\n  Hermes-3 face value / SOLPS guard (= target) value:")
    face_matches_solps_guard = True
    for variable_name in ["Ne", "Te", "Td+", "Vd+"]:
        face_value_hermes  = hermes_face_value(hermes_last_cell, hermes_guard_cell, variable_name)
        guard_value_solps  = solps_value(variable_name, solps_guard_index)
        ratio_face_to_guard = np.abs(face_value_hermes) / np.abs(guard_value_solps)
        print_ratio_summary(f"{variable_name} face(H3)/guard(SOLPS)", ratio_face_to_guard)
        if not np.all(np.abs(ratio_face_to_guard - 1) < tolerance):
            face_matches_solps_guard = False

    # -----------------------------------------------------------------------
    # E. Sheath-extrapolated face values (sheath_boundary_simple, mode 1)
    #    Commented out: only meaningful when sheath_boundary_simple is in
    #    [hermes] components. The frozen-plasma runs do not use it.
    # -----------------------------------------------------------------------
    # print_section_header("E. Hermes-3 sheath-extrapolated face value vs SOLPS", "-")
    #
    # density_face_extrapolated = sheath_extrapolated_face_value(
    #     hermes_last_cell, hermes_second_last_cell, "Nd+")
    # electron_temperature_face_extrapolated = sheath_extrapolated_face_value(
    #     hermes_last_cell, hermes_second_last_cell, "Te")
    # ion_temperature_face_extrapolated = sheath_extrapolated_face_value(
    #     hermes_last_cell, hermes_second_last_cell, "Td+")
    # sound_speed_face = np.sqrt(ELECTRON_CHARGE
    #                            * (ion_temperature_face_extrapolated + electron_temperature_face_extrapolated)
    #                            / (ion_atomic_mass * PROTON_MASS))
    # # sheath boundary forces the face velocity to be at least the sound speed
    # velocity_face_extrapolated = np.maximum(np.abs(hermes_last_cell["Vd+"].values), sound_speed_face)
    #
    # print_ratio_summary("n  solver/SOLPS",
    #                     density_face_extrapolated / solps_value("Ne", solps_guard_index))
    # print_ratio_summary("Te solver/SOLPS",
    #                     electron_temperature_face_extrapolated / solps_value("Te", solps_guard_index))
    # print_ratio_summary("Ti solver/SOLPS",
    #                     ion_temperature_face_extrapolated / solps_value("Td+", solps_guard_index))
    #
    # flow_rebuilt_from_extrapolation = float(
    #     (density_face_extrapolated * velocity_face_extrapolated * face_area_hermes).sum())
    # target_cells = ds.hermesm.select_region(f"{target_name}_target").squeeze()
    # flow_from_recycle_source = (float((target_cells["Sd_target_recycle"] * target_cells["dv"]).sum())
    #                             / recycle_multiplier)
    # print(f"\n  Particle flow [s^-1]:  rebuilt n*cs*A {flow_rebuilt_from_extrapolation:.4e}   "
    #       f"from Sd_target_recycle {flow_from_recycle_source:.4e}   "
    #       f"SOLPS {bal.target_flux[target_name]:.4e}")

    # -----------------------------------------------------------------------
    # F. Verdict: compare the TARGET value only
    #    SOLPS target = guard cell (thin, see C)
    #    Hermes-3 target = face = (last + guard)/2 (thick guard, see D)
    # -----------------------------------------------------------------------
    print_section_header("F. Verdict: is the target value sonic in each code?", "-")
    sonic_rings_solps_target,  total_rings = count_sonic_rings(mach_solps_guard)
    sonic_rings_hermes_target, _           = count_sonic_rings(mach_hermes_face)

    # np.ma.median: the SOLPS arrays are numpy masked arrays, plain np.median warns about them
    verdict_table = pd.DataFrame(
        [["SOLPS (target = guard)",   sonic_rings_solps_target,  np.ma.median(mach_solps_guard)],
         ["Hermes-3 (target = face)", sonic_rings_hermes_target, np.ma.median(mach_hermes_face)]],
        columns=["code", f"sonic target /{total_rings}", "median M target"],
    ).set_index("code")
    print(verdict_table.to_string(float_format="{:.3f}".format))
    print(f"  (sonic = |M - 1| < {tolerance})\n")

    solps_target_is_sonic  = sonic_rings_solps_target  > 0.9 * total_rings
    hermes_target_is_sonic = sonic_rings_hermes_target > 0.9 * total_rings

    verdict_lines = []
    if solps_target_is_sonic:
        verdict_lines.append("✅ SOLPS: guard cell is sonic → guard value = target value")
    else:
        verdict_lines.append("❌ SOLPS: guard cell is NOT sonic (check table C)")

    if hermes_target_is_sonic:
        verdict_lines.append("✅ Hermes-3: face (last+guard)/2 is sonic → face value = target value")
    else:
        verdict_lines.append("❌ Hermes-3: face is NOT sonic (check table D)")

    if face_matches_solps_guard:
        verdict_lines.append(f"✅ Hermes-3 face matches SOLPS guard "
                             f"(Ne, Te, Td+, Vd+ within ±{tolerance:.0%})")
    else:
        verdict_lines.append("❌ Hermes-3 face does NOT match SOLPS guard (check ratios in D)")

    for line in verdict_lines:
        print("  " + line)
    verdicts_per_target[target_name] = verdict_lines


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print_section_header("SUMMARY")
print(f"Hermes-3 case : {args.input}")
print(f"SOLPS case    : {args.solps}")
print(f"Components    : {ds.options['hermes']['components']}")
for target_name, verdict_lines in verdicts_per_target.items():
    print(f"\n[{target_name}]")
    for line in verdict_lines:
        print("  " + line)

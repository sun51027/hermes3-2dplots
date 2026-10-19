"""
Quick health check of one (or more) Hermes-3 cases, with a Markdown + PDF report.

Script version of calculator.ipynb. It reads BOTH the input file (BOUT.inp)
and the result file (BOUT.squash.nc / BOUT.dmp.*.nc), runs a list of checks,
and writes for every case:

    check_reports/<case_id>/<case_id>_check_report.md    (open in VS Code preview)
    check_reports/<case_id>/<case_id>_check_report.pdf   (made from the .md by xelatex)
    check_reports/<case_id>/figures/fig_XX.png
    check_reports/<case_id>/<case_id>_report_history.csv  (version 1, 2, 3, ... with dates)

This file holds the settings and the checks. Printing, Markdown, PDF and the version
history are in check_hermes_case_output.py (same folder, normally no need to edit it).

Every section of the report has: the formulas, the variables used, and result tables.
Every PDF page has a footer: case name, page number, report version and date.
Delete the history CSV to restart the version count at 1.

Usage (from hermes3-2dplots/):
    python check_hermes_case.py -i 261001_highrsl_fixplsbg_fromsolps_cx100_1e21_roomT_gd
    python check_hermes_case.py -i caseA caseB            # one report per case
    python check_hermes_case.py -i caseA --no-plots       # tables only (fast)
    python check_hermes_case.py -i caseA --expected-puff-rate 1e21

Use from another script or notebook:
    from check_hermes_case import load_case, run_checks, write_report
    case = load_case("261001_highrsl_fixplsbg_fromsolps_cx100_1e21_roomT_gd")
    results = run_checks(case)
    write_report(case, results)

    # or call one check / helper alone:
    from check_hermes_case import integrate_over_region, get_input_number
    print(integrate_over_region(case.ds, "Sd_src", "all_noguards"))
    print(get_input_number(case.input_options, "Pe:pressure_source"))

------------------------------------------------------------------------------
HOW TO ADD THINGS (no new code needed for 1-4)
------------------------------------------------------------------------------
1. Show one more value from BOUT.inp
       -> add "section:key" to INPUT_VALUES_TO_SHOW           (e.g. "d:neutral_lmax")
2. Show one more volume integral sum(variable * dv) from the result file
       -> add one dict to VOLUME_INTEGRALS_TO_SHOW
3. Compare a number from BOUT.inp with a number from the result file
       -> add one dict to INPUT_VS_RESULT_COMPARISONS
4. Add a 2D polygon plot  -> add the variable to POLYGON_PAGES
   Add a 1D profile line  -> add the variable to PROFILE_ROWS
5. A completely new check:
       a) copy the template "check_template_example" near the end of section 5
       b) fill it in: add_formula(), add_variable(), add_table(), add_figure()
       c) add the function name to CHECKS_TO_RUN (section 7)

Status words used in the report:
    PASS = calculated value matches the expected value within PASS_TOLERANCE
    WARN = within WARN_TOLERANCE
    FAIL = worse than WARN_TOLERANCE, or the check crashed
    INFO = nothing to compare, information only
"""

import argparse
import re
import sys
import traceback
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

sys.path.append(r"/Users/zero/workspace/phd_work/hms_output/sdtools")
import xbout                     # noqa: F401  registers the ds["var"].bout.polygon() plotting
import hermes3.accessors         # noqa: F401  registers ds.hermesm.select_region()
from boutdata.data import BoutOptionsFile
from hermes3.accessors import _select_region   # returns the (x, theta) indices of a named region
from hermes3.case_db import CaseDB
from hermes3.selectors import get_1d_poloidal_data, get_1d_radial_data
from hermes3.utils import constants

# Output part (terminal, Markdown, PDF, version history) lives in check_hermes_case_output.py,
# in the same folder as this script.
sys.path.append(str(Path(__file__).resolve().parent))
from check_hermes_case_output import (STATUS_RANK, format_number, format_status,
                                      print_results, write_report_files)


# ============================================================================
# 1. USER SETTINGS  (edit these freely)
# ============================================================================
DEFAULT_CASE_ID  = "261001_highrsl_fixplsbg_fromsolps_cx100_1e21_roomT_gd"
CASE_DIRECTORY   = r"/Users/zero/workspace/phd_work/hms_output"
GRID_DIRECTORY   = r"/Users/zero/workspace/phd_work/hms_output/gridfile"
REPORT_DIRECTORY = Path(__file__).resolve().parent / "check_reports"

PASS_TOLERANCE = 0.01   # relative difference below this -> PASS  (1 %)
WARN_TOLERANCE = 0.10   # relative difference below this -> WARN  (10 %), otherwise FAIL

ELECTRON_CHARGE = constants("q_e")   # [C] = [J/eV]

# Sections in BOUT.inp that hold a core power source:
#   power_in [MW], core_volume [m^3], source = H(x) * H(b - x) * pressure_source
# The matching result variable is "<section>_src", e.g. "Pe_src".
CORE_SOURCE_SECTIONS = ["Pe", "Pd+"]

# Gas puff: density and pressure source variables in the result file
PUFF_DENSITY_SOURCE  = "Sd_src"     # [1/s/m^3]
PUFF_PRESSURE_SOURCE = "Pd_src"     # [Pa/s], energy source = 3/2 * Pd_src

# Optional expected puff values. None = do not compare, just show.
EXPECTED_PUFF_RATE       = None     # [atoms/s], e.g. 1e21
EXPECTED_PUFF_ENERGY_EV  = None     # [eV per atom], e.g. 0.02533 (room temperature)

# Section in BOUT.inp that holds target_recycle_multiplier
RECYCLING_ION_SECTION = "d+"

# Core ionising (neutral_boundary with d:ionising_core = true): neutrals that reach
# the core boundary are removed (neutral sink) and put back as ions (ion source).
CORE_IONISING_NEUTRAL_SINK = "Nd_core_sink"     # [1/s/m^3], negative
CORE_IONISING_ION_SOURCE   = "Nd+_core_source"  # [1/s/m^3], positive


# ---- 1a. Values from BOUT.inp to list in the report --------------------------
# Format "section:key". Use just "key" for values at the top of BOUT.inp.
INPUT_VALUES_TO_SHOW = [
    "nout",
    "timestep",
    "mesh:file",
    "hermes:components",
    # "hermes:Nnorm",
    "d+:type",
    "e:type",
    "d:type",
    "Pe:power_in",
    "Pe:core_volume",
    "Pe:pressure_source",
    "Pe:source",
    "Pd+:power_in",
    "Pd+:core_volume",
    "Pd+:pressure_source",
    "Pd+:source",
    # "d+:target_recycle_multiplier",
    # "d+:target_recycle_energy",
    "d+:neutral_pump",
    # "d+:pump_recycle_multiplier",
    "d:K_cx_multiplier",
    "d:flux_limit",
    "d:neutral_lmax",
    # "d:ionising_core",
    # "d:ionising_core_iz_energy_loss",
]


# ---- 1b. Volume integrals from the result file -----------------------------
# Each line computes   factor * sum(variable * dv)   over a named region.
# Region names: "all_noguards", "core_noguards", "core_edge" (= first core ring),
#               "sol", "outer_lower_target", ...  (see sdtools/hermes3/accessors.py)
# variable = None  -> just the volume of the region.
# Units may use LaTeX math, e.g. "m$^{-3}$".
VOLUME_INTEGRALS_TO_SHOW = [
    {"label": "Neutral atoms in whole domain", "variable": "Nd",     "region": "all_noguards",  "factor": 1.0, "unit": "atoms"},
    {"label": "Neutral atoms in core",         "variable": "Nd",     "region": "core_noguards", "factor": 1.0, "unit": "atoms"},
    {"label": "Ions in whole domain",          "variable": "Nd+",    "region": "all_noguards",  "factor": 1.0, "unit": "ions"},
    {"label": "Ionisation in core",            "variable": "Sd+_iz", "region": "core_noguards", "factor": 1.0, "unit": "s$^{-1}$"},
    {"label": "CX energy transfer",            "variable": "Edd+_cx","region": "all_noguards",  "factor": 1.0, "unit": "W"},
]


# ---- 1c. Compare a BOUT.inp number with a result-file number ---------------
# expected   = input_factor  * (value of input_key in BOUT.inp)
# calculated = result_factor * sum(result_variable * dv) over result_region
INPUT_VS_RESULT_COMPARISONS = [
    {"label":           "[Pe] core_volume vs volume of first core ring",
     "input_key":       "Pe:core_volume", "input_factor":  1.0,
     "result_variable": None,             "result_region": "core_edge", "result_factor": 1.0,
     "unit":            "m$^3$"},
    # Template: remove the "#" to switch it on
    # {"label":           "[Pd+] power_in vs 3/2 * sum(Pd+_src * dv)",
    #  "input_key":       "Pd+:power_in",  "input_factor":  1e6,
    #  "result_variable": "Pd+_src",       "result_region": "all_noguards", "result_factor": 1.5,
    #  "unit":            "W"},
]


# ---- 1d. 2D polygon plots ---------------------------------------------------
# One dict = one figure. Options per variable:
#   "logscale": True/False, "absolute": True -> plot |value|
POLYGON_PAGES = [
    {"page_title": "Sources, sinks and pump",
     "variables": [
         {"name": "Sd_src",            "logscale": True},
         {"name": "Pd_src",            "logscale": True},
         {"name": "Sd_pump"},
         {"name": "is_pump"},
         {"name": "Sd_target_recycle", "logscale": True, "absolute": True},
         {"name": "Pe_src"},
         {"name": "Pd+_src"},
         {"name": "Sd+_iz",            "logscale": True, "absolute": True},
     ]},
    {"page_title": "Anomalous transport coefficients",
     "variables": [
         {"name": "anomalous_Chi_d+", "logscale": True},
         {"name": "anomalous_Chi_e"},
         {"name": "anomalous_D_d+"},
         {"name": "anomalous_D_e"},
         {"name": "anomalous_nu_d+"},
         {"name": "anomalous_nu_e"},
     ]},
    {"page_title": "Neutral state variables",
     "variables": [
         {"name": "Td", "absolute": True},
         {"name": "Nd", "absolute": True, "logscale": True},
         {"name": "Pd", "absolute": True, "logscale": True},
         {"name": "Vd"},
     ]},
    {"page_title": "Plasma state variables",
     "variables": [
         {"name": "Te"},
         {"name": "Td+"},
         {"name": "Ne",  "logscale": True},
         {"name": "Vd+"},
     ]},
    {"page_title": "Where is the solver struggling? |ddt(...)|",
     "variables": [
         {"name": "ddt(Nd)",   "absolute": True, "logscale": True},
         {"name": "ddt(Pd)",   "absolute": True, "logscale": True},
         {"name": "ddt(NVd)",  "absolute": True, "logscale": True},
         {"name": "ddt(Nd+)",  "absolute": True, "logscale": True},
         {"name": "ddt(Pd+)",  "absolute": True, "logscale": True},
         {"name": "ddt(Pe)",   "absolute": True, "logscale": True},
     ]},
]


# ---- 1e. 1D profiles (one figure: rows x columns) ---------------------------
# Rows = what to plot. Columns = where to take the profile.
PROFILE_ROWS = [
    {"ylabel": "Temperature [eV]", "variables": ["Te", "Td+", "Td"], "yscale": "linear"},
    {"ylabel": "Density [m^-3]",   "variables": ["Ne", "Nd+", "Nd"], "yscale": "log"},
    {"ylabel": "Pressure [Pa]",    "variables": ["Pe", "Pd+", "Pd"], "yscale": "log"},
]
# kind = "poloidal": along one SOL ring. sepadd = 0 is the FIRST SOL ring
#                    (sepadd = 1 is the second one).
# kind = "radial"  : across the flux surfaces at region "omp" or "imp".
PROFILE_COLUMNS = [
    {"title": "First SOL ring, outer lower leg", "kind": "poloidal",
     "region": "outer_lower_sol", "sepadd": 0, "target_first": True},
    {"title": "Outer midplane", "kind": "radial", "region": "omp"},
]


# ============================================================================
# 2. DATA CONTAINERS
# ============================================================================
class CaseData:
    """
    Everything the checks need, in one place (like a C++ struct).

        case.case_id         name of the case folder
        case.ds_all_times    xarray Dataset with all time steps
        case.ds              xarray Dataset at the last time step (same as the notebook)
        case.input_options   BOUT.inp read by boutdata (sections and keys)
        case.input_file_path path of BOUT.inp
        case.result_file_path path of BOUT.squash.nc (or BOUT.dmp.0.nc)
    """

    def __init__(self, case_id, ds_all_times, input_options, input_file_path, result_file_path):
        self.case_id = case_id
        self.ds_all_times = ds_all_times
        self.ds = ds_all_times.isel(t=-1)
        self.input_options = input_options
        self.input_file_path = input_file_path
        self.result_file_path = result_file_path


# Where a variable comes from (used in the "Variables used" tables)
FROM_RESULT = "result file (nc)"
FROM_INPUT  = "BOUT.inp"
FROM_GRID   = "grid metadata"


class CheckResult:
    """
    Output of one check (one section of the report).

        title, description   section heading and one-line explanation
        status               INFO / PASS / WARN / FAIL (only gets worse, never better)
        formulas             LaTeX formulas, shown at the top of the section
        variables            rows of the "Variables used" table
        blocks               results in order: text, tables and figures

    Text may use Markdown: `code`, **bold**, $math$.
    """

    def __init__(self, title, description=""):
        self.title = title
        self.description = description
        self.status = "INFO"
        self.formulas = []
        self.variables = []
        self.blocks = []

    def update_status(self, new_status):
        if STATUS_RANK[new_status] > STATUS_RANK[self.status]:
            self.status = new_status

    def add_formula(self, latex_formula):
        """One display formula, LaTeX syntax without the $$."""
        self.formulas.append(latex_formula)

    def add_variable(self, name, source, meaning, unit=""):
        """One row of the 'Variables used' table."""
        self.variables.append({"Variable": f"`{name}`", "From": source,
                               "Meaning": meaning, "Unit": unit})

    def add_text(self, text):
        self.blocks.append({"type": "text", "text": text})

    def add_table(self, rows, caption=""):
        """rows = list of dicts, one dict per row, same keys in every dict (= column names)."""
        self.blocks.append({"type": "table", "rows": rows, "caption": caption})

    def add_figure(self, figure, caption=""):
        self.blocks.append({"type": "figure", "figure": figure, "caption": caption})


# ============================================================================
# 3. LOADING, READING AND FORMATTING HELPERS
# ============================================================================
def load_case(case_id, case_directory=CASE_DIRECTORY, grid_directory=GRID_DIRECTORY,
              use_squash=True):
    """Load the result file and BOUT.inp of one case. Returns a CaseData."""
    db = CaseDB(case_dir=case_directory, grid_dir=grid_directory)
    case = db.load_case_2D(case_id, use_squash=use_squash, verbose=True)
    case.extract_2d_tokamak_geometry()

    case_folder = Path(db.casepaths[case_id])
    input_file_path = case_folder / "BOUT.inp"
    input_options = read_input_file(input_file_path)

    squash_file_path = case_folder / "BOUT.squash.nc"
    if use_squash and squash_file_path.exists():
        result_file_path = squash_file_path
    else:
        result_file_path = case_folder / "BOUT.dmp.0.nc"
    return CaseData(case_id, case.ds, input_options, input_file_path, result_file_path)


def read_input_file(input_file_path):
    """Read BOUT.inp. The warning about x/y/z arrays is harmless and hidden."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        input_options = BoutOptionsFile(str(input_file_path))
    return input_options


def split_input_key(full_key):
    """'Pe:power_in' -> ('Pe', 'power_in').  'nout' -> ('', 'nout')."""
    if ":" in full_key:
        section_name, key_name = full_key.rsplit(":", 1)
    else:
        section_name = ""
        key_name = full_key
    return section_name, key_name


def get_input_section(input_options, section_name):
    """Return one section of BOUT.inp, or None if it does not exist."""
    if section_name == "":
        return input_options
    try:
        return input_options[section_name]
    except KeyError:
        return None


def get_input_raw(input_options, full_key):
    """Value exactly as written in BOUT.inp (number or text). None if not set."""
    section_name, key_name = split_input_key(full_key)
    section = get_input_section(input_options, section_name)
    if section is None:
        return None
    try:
        return section[key_name]
    except KeyError:
        return None


def get_input_number(input_options, full_key):
    """
    Value from BOUT.inp as a float. Expressions such as
        pressure_source = power_in * 1e6 * 2/3 / core_volume
    are evaluated. Returns None if not set or not a plain number (e.g. contains H(x)).
    """
    raw_value = get_input_raw(input_options, full_key)
    if raw_value is None:
        return None
    if isinstance(raw_value, (int, float)):
        return float(raw_value)

    section_name, key_name = split_input_key(full_key)
    section = get_input_section(input_options, section_name)
    try:
        return float(section.evaluate_scalar(key_name))
    except Exception:
        return None


def get_input_bool(input_options, full_key, default_value=False):
    """Value from BOUT.inp as True/False."""
    raw_value = get_input_raw(input_options, full_key)
    if raw_value is None:
        return default_value
    return str(raw_value).strip().lower() in ["true", "yes", "y", "1"]


def integrate_over_region(ds, variable_name, region_name="all_noguards", factor=1.0):
    """
    factor * sum(variable * dv) over a named region, at one time step.
    variable_name = None -> volume of the region.
    Returns None if the variable is not in the result file.
    """
    if variable_name is not None and variable_name not in ds:
        return None
    region = ds.hermesm.select_region(region_name)
    if variable_name is None:
        integrand = region["dv"]
    else:
        integrand = region[variable_name] * region["dv"]
    return factor * float(integrand.sum())


def sum_over_domain_versus_time(ds_all_times, variable_name):
    """sum(variable * dv) over the domain without guards, for every time step."""
    if variable_name not in ds_all_times:
        return None
    domain = ds_all_times.hermesm.select_region("all_noguards")
    return (domain[variable_name] * domain["dv"]).sum(["x", "theta"])


def find_heaviside_window(source_expression):
    """
    Read the x range of a source like  H(x) * H(0.035 - x) * pressure_source.
    Returns (lower_bound, upper_bound). A missing bound is -inf or +inf.
    Returns None if the expression is not text.
    """
    if not isinstance(source_expression, str):
        return None
    expression = source_expression.replace(" ", "").lower()
    lower_bound = -np.inf
    upper_bound = np.inf

    if "h(x)" in expression:
        lower_bound = 0.0
    lower_match = re.search(r"h\(x-([0-9.e+\-]+)\)", expression)
    if lower_match is not None:
        lower_bound = float(lower_match.group(1))
    upper_match = re.search(r"h\(([0-9.e+\-]+)-x\)", expression)
    if upper_match is not None:
        upper_bound = float(upper_match.group(1))
    return lower_bound, upper_bound


def judge_difference(calculated_value, expected_value):
    """Returns (status, relative_difference). Status INFO if a value is missing."""
    if calculated_value is None or expected_value is None:
        return "INFO", None
    if expected_value == 0:
        relative_difference = abs(calculated_value)
    else:
        relative_difference = abs(calculated_value - expected_value) / abs(expected_value)

    if relative_difference <= PASS_TOLERANCE:
        status = "PASS"
    elif relative_difference <= WARN_TOLERANCE:
        status = "WARN"
    else:
        status = "FAIL"
    return status, relative_difference


def make_comparison_row(result, label, calculated_value, expected_value, unit=""):
    """
    One table row comparing a calculated value with an expected value.
    Also updates the status of the whole check.
    """
    status, relative_difference = judge_difference(calculated_value, expected_value)
    result.update_status(status)
    if relative_difference is None:
        difference_text = "n/a"
    else:
        difference_text = f"{relative_difference:.2%}"
    return {"Item": label,
            "Calculated": format_number(calculated_value),
            "Expected": format_number(expected_value),
            "Unit": unit,
            "Difference": difference_text,
            "Status": format_status(status)}


# ============================================================================
# 4. GRID / RING HELPERS
# ============================================================================
def get_radial_grid_sizes(ds):
    """
    Returns (number_of_x_guards, separatrix_index, number_of_radial_cells).
        separatrix_index       = ixseps1 = x index of the first SOL ring (guards counted)
        number_of_radial_cells = radial cells without guards (BOUT++ calls this MX)
    """
    number_of_x_guards = int(ds.metadata["MXG"])
    separatrix_index = int(ds.metadata["ixseps1"])
    number_of_radial_cells = ds.sizes["x"] - 2 * number_of_x_guards
    return number_of_x_guards, separatrix_index, number_of_radial_cells


def build_radial_ring_table(ds):
    """
    One row per radial ring (x index), with the BOUT++ "x" used in BOUT.inp expressions.

    BOUT++ (boutmesh.cxx, GlobalX with symmetricGlobalX = true):
        x = (i + 0.5) / MX        i = 0 for the first cell next to the core boundary
    So x depends on the cell INDEX only, not on the cell width dx.
    "old notebook x" = cumsum(dx)/sum(dx) at the OMP, as in calculator.ipynb, for comparison.
    """
    number_of_x_guards, separatrix_index, number_of_radial_cells = get_radial_grid_sizes(ds)
    total_x_points = ds.sizes["x"]

    core_theta_indices = _select_region(ds, "core_noguards")[1]
    domain_theta_indices = _select_region(ds, "all_noguards")[1]
    cell_volume = ds["dv"].transpose("x", "theta")
    core_ring_volumes = cell_volume.isel(theta=core_theta_indices).sum("theta").values
    full_ring_volumes = cell_volume.isel(theta=domain_theta_indices).sum("theta").values

    # index of omp_profile = x index with guards counted (guards themselves removed)
    omp_profile = get_1d_radial_data(ds, ["dx"], region="omp")
    old_notebook_x = omp_profile["dx"].cumsum() / omp_profile["dx"].sum()

    table_rows = []
    for x_index in range(number_of_x_guards, total_x_points - number_of_x_guards):
        domain_index = x_index - number_of_x_guards
        if x_index < separatrix_index:
            ring_label = f"core {domain_index + 1}"
            is_core_ring = True
            ring_volume = core_ring_volumes[x_index]
        else:
            ring_label = f"SOL {x_index - separatrix_index + 1}"
            is_core_ring = False
            ring_volume = full_ring_volumes[x_index]

        table_rows.append({
            "x index": x_index,
            "ring": ring_label,
            "is core": is_core_ring,
            "x inner face": domain_index / number_of_radial_cells,
            "x centre": (domain_index + 0.5) / number_of_radial_cells,
            "x outer face": (domain_index + 1) / number_of_radial_cells,
            "old notebook x": old_notebook_x.get(x_index, np.nan),
            "Srad OMP [mm]": 1000 * omp_profile["Srad"].get(x_index, np.nan),
            "ring volume": ring_volume,
        })
    return pd.DataFrame(table_rows).set_index("x index")


def select_rings_to_show(ring_table, separatrix_index):
    """First 5 rings, and 2 rings on each side of the separatrix."""
    wanted_indices = list(ring_table.index[:5])
    for offset in [-2, -1, 0, 1, 2]:
        wanted_indices.append(separatrix_index + offset)
    rows_to_keep = []
    for x_index in ring_table.index:
        if x_index in wanted_indices:
            rows_to_keep.append(x_index)
    return ring_table.loc[rows_to_keep]


def find_rings_in_window(ring_table, lower_bound, upper_bound, core_only):
    """x indices of rings whose centre is inside lower_bound < x < upper_bound."""
    inside_window = (ring_table["x centre"] > lower_bound) & (ring_table["x centre"] < upper_bound)
    if core_only:
        inside_window = inside_window & ring_table["is core"]
    return list(ring_table.index[inside_window])


def find_rings_with_positive_values(ds, variable_name):
    """x indices where the variable is > 0 in at least one poloidal cell."""
    values = ds[variable_name].transpose("x", "theta").values
    ring_has_positive_value = (values > 0).any(axis=1)
    return list(np.where(ring_has_positive_value)[0])


def ring_labels_as_text(ring_table, x_indices):
    """[37, 38] -> 'core 1, core 2'. Empty list -> 'none'."""
    if len(x_indices) == 0:
        return "none"
    labels = []
    for x_index in x_indices:
        if x_index in ring_table.index:
            labels.append(ring_table.loc[x_index, "ring"])
        else:
            labels.append(f"guard x={x_index}")
    return ", ".join(labels)


# ============================================================================
# 5. CHECKS  (each returns one CheckResult = one section of the report)
# ============================================================================
def check_input_file_values(case):
    """List the BOUT.inp values in INPUT_VALUES_TO_SHOW."""
    result = CheckResult("A. Input file values (BOUT.inp)",
                         "Values as written in BOUT.inp. 'Evaluated' = expression with the "
                         "other keys substituted (boutdata `evaluate_scalar`).")
    result.add_formula(r"p_{\rm src} = \frac{2}{3}\,\frac{P_{\rm in}\,[\mathrm{MW}] \times 10^{6}}{V_{\rm core}}"
                       r"\quad [\mathrm{Pa\,s^{-1}}]")
    result.add_formula(r"\text{source} = H(x - a)\,H(b - x)\,p_{\rm src}")
    result.add_variable("<section>:<key>", FROM_INPUT, "every key in INPUT_VALUES_TO_SHOW")

    table_rows = []
    for full_key in INPUT_VALUES_TO_SHOW:
        raw_value = get_input_raw(case.input_options, full_key)
        evaluated_value = get_input_number(case.input_options, full_key)
        if raw_value is None:
            raw_text = "(not set)"
        else:
            raw_text = f"`{raw_value}`"
        if evaluated_value is None or isinstance(raw_value, (int, float)):
            evaluated_text = ""
        else:
            evaluated_text = format_number(evaluated_value)
        table_rows.append({"Key": f"`{full_key}`", "Value in BOUT.inp": raw_text,
                           "Evaluated": evaluated_text})
    result.add_table(table_rows)
    return result


def check_total_power_source(case):
    """
    Core power from the result file:  P = 3/2 * sum(P_src * dv)   [W]
    compared with power_in [MW] in BOUT.inp, for every section in CORE_SOURCE_SECTIONS.
    """
    result = CheckResult("B. Total power source",
                         "Power put into the core by the pressure sources, compared with "
                         "`power_in` in BOUT.inp. Sum over the domain without guard cells.")
    result.add_formula(r"P_s = \frac{3}{2} \sum_{i \in \text{domain}} S_{P,s,i}\,\Delta V_i \quad [\mathrm{W}]")
    result.add_formula(r"P_s^{\rm expected} = P_{{\rm in},s} \times 10^{6}")
    result.add_formula(r"P_{\rm puff} = \frac{3}{2} \sum_{i} S_{P,d,i}\,\Delta V_i")
    for section_name in CORE_SOURCE_SECTIONS:
        result.add_variable(f"{section_name}_src", FROM_RESULT,
                            f"pressure source $S_P$ of {section_name}", "Pa s$^{-1}$")
        result.add_variable(f"{section_name}:power_in", FROM_INPUT,
                            f"input power $P_{{\\rm in}}$ of {section_name}", "MW")
    result.add_variable(PUFF_PRESSURE_SOURCE, FROM_RESULT, "neutral pressure source (puff)", "Pa s$^{-1}$")
    result.add_variable("dv", FROM_RESULT, r"cell volume $\Delta V$", "m$^3$")

    ds = case.ds
    table_rows = []
    total_calculated_power = 0.0
    total_expected_power = 0.0
    number_of_sources_found = 0

    for section_name in CORE_SOURCE_SECTIONS:
        source_variable = f"{section_name}_src"
        power_in_megawatt = get_input_number(case.input_options, f"{section_name}:power_in")
        calculated_power = integrate_over_region(ds, source_variable, "all_noguards", factor=1.5)

        if power_in_megawatt is None:
            expected_power = None
        else:
            expected_power = power_in_megawatt * 1e6

        if calculated_power is None:
            result.add_text(f"`{source_variable}` is not in the result file "
                            f"(frozen plasma or no source): no comparison for [{section_name}].")
        else:
            number_of_sources_found += 1
            total_calculated_power += calculated_power
            if expected_power is not None:
                total_expected_power += expected_power
        table_rows.append(make_comparison_row(result, f"[{section_name}] power", calculated_power,
                                              expected_power, "W"))

    if number_of_sources_found > 0:
        table_rows.append(make_comparison_row(result, "**Total core power**", total_calculated_power,
                                              total_expected_power, "W"))

    neutral_heating = integrate_over_region(ds, PUFF_PRESSURE_SOURCE, "all_noguards", factor=1.5)
    table_rows.append(make_comparison_row(result, "Neutral heating by puff $P_{\\rm puff}$",
                                          neutral_heating, None, "W"))
    result.add_table(table_rows)
    return result


def check_gas_puff_source(case):
    """Where the puff is, the total puff rate and the energy per injected atom."""
    result = CheckResult("C. Gas puff source",
                         "Cells where the puff density source is > 0, the total puff rate and "
                         "the energy given to each injected atom.")
    result.add_formula(r"\Gamma_{\rm puff} = \sum_i S_{n,i}\,\Delta V_i \quad [\mathrm{s^{-1}}]")
    result.add_formula(r"P_{\rm puff} = \frac{3}{2}\sum_i S_{P,i}\,\Delta V_i \quad [\mathrm{W}]")
    result.add_formula(r"E_{\rm atom} = \frac{P_{\rm puff}}{e\,\Gamma_{\rm puff}} \quad [\mathrm{eV}]"
                       r"\qquad (= T \text{ in calculator.ipynb})")
    result.add_variable(PUFF_DENSITY_SOURCE, FROM_RESULT, "puff density source $S_n$", "m$^{-3}$ s$^{-1}$")
    result.add_variable(PUFF_PRESSURE_SOURCE, FROM_RESULT, "puff pressure source $S_P$", "Pa s$^{-1}$")
    result.add_variable("dv", FROM_RESULT, r"cell volume $\Delta V$", "m$^3$")
    result.add_variable("R, Z", FROM_RESULT, "cell centre position", "m")
    result.add_variable("q_e", "constant", "electron charge $e$", "C")

    ds = case.ds
    if PUFF_DENSITY_SOURCE not in ds:
        result.add_text(f"`{PUFF_DENSITY_SOURCE}` is not in the result file.")
        return result

    density_source = ds[PUFF_DENSITY_SOURCE].transpose("x", "theta")
    x_indices, theta_indices = np.where(density_source.values > 0)
    if len(x_indices) == 0:
        result.add_text(f"**WARN** `{PUFF_DENSITY_SOURCE}` is zero everywhere: no puff.")
        result.update_status("WARN")
        return result

    has_pressure_source = PUFF_PRESSURE_SOURCE in ds
    cell_rows = []
    for x_index, theta_index in zip(x_indices, theta_indices):
        one_cell = ds.isel(x=x_index, theta=theta_index)
        row = {"x index": int(x_index), "theta index": int(theta_index),
               "R [m]": format_number(float(one_cell["R"])),
               "Z [m]": format_number(float(one_cell["Z"])),
               "$S_n$ [m$^{-3}$ s$^{-1}$]": format_number(float(one_cell[PUFF_DENSITY_SOURCE])),
               r"$\Delta V$ [m$^3$]": format_number(float(one_cell["dv"]))}
        if has_pressure_source:
            row["$S_P$ [Pa s$^{-1}$]"] = format_number(float(one_cell[PUFF_PRESSURE_SOURCE]))
        cell_rows.append(row)
    result.add_table(cell_rows, f"Puff cells ({len(cell_rows)})")

    # same sum as the notebook: over all cells where the source is non-zero
    puff_rate = float((ds[PUFF_DENSITY_SOURCE] * ds["dv"]).sum())
    summary_rows = [make_comparison_row(result, r"Puff rate $\Gamma_{\rm puff}$", puff_rate,
                                        EXPECTED_PUFF_RATE, "s$^{-1}$")]
    if has_pressure_source:
        puff_heating_power = 1.5 * float((ds[PUFF_PRESSURE_SOURCE] * ds["dv"]).sum())
        energy_per_atom_ev = puff_heating_power / (puff_rate * ELECTRON_CHARGE)
        summary_rows.append(make_comparison_row(result, r"Puff heating $P_{\rm puff}$",
                                                puff_heating_power, None, "W"))
        summary_rows.append(make_comparison_row(result, r"Energy per atom $E_{\rm atom}$",
                                                energy_per_atom_ev, EXPECTED_PUFF_ENERGY_EV, "eV"))
    result.add_table(summary_rows, "Puff totals")
    return result


def check_plasma_volumes(case):
    """Volume of the whole domain, core, SOL, PFR and the first core ring."""
    result = CheckResult("D. Plasma volumes", "Sum of cell volumes per region, guard cells excluded.")
    result.add_formula(r"V_{\rm region} = \sum_{i \in \text{region}} \Delta V_i")
    result.add_formula(r"V_{\rm SOL} = V_{\rm domain} - V_{x < x_{\rm sep}}, \qquad "
                       r"V_{\rm PFR} = V_{x < x_{\rm sep}} - V_{\rm core}")
    result.add_variable("dv", FROM_RESULT, r"cell volume $\Delta V$", "m$^3$")
    result.add_variable("MXG", FROM_GRID, "number of x guard cells", "")
    result.add_variable("ixseps1", FROM_GRID, "x index of the first SOL ring", "")

    ds = case.ds
    number_of_x_guards, separatrix_index, _ = get_radial_grid_sizes(ds)

    domain_volume = integrate_over_region(ds, None, "all_noguards")
    core_volume = integrate_over_region(ds, None, "core_noguards")
    first_core_ring_volume = integrate_over_region(ds, None, "core_edge")

    # inside the separatrix = core + PFR
    domain = ds.hermesm.select_region("all_noguards")
    inside_separatrix_volume = float(
        domain["dv"].isel(x=slice(0, separatrix_index - number_of_x_guards)).sum())
    sol_volume = domain_volume - inside_separatrix_volume
    pfr_volume = inside_separatrix_volume - core_volume

    table_rows = [
        {"Region": "Whole domain",    "How it is selected": "`all_noguards`",
         "Volume [m$^3$]": format_number(domain_volume)},
        {"Region": "Core",            "How it is selected": "`core_noguards`",
         "Volume [m$^3$]": format_number(core_volume)},
        {"Region": "SOL",             "How it is selected": "domain cells with x index $\\ge$ ixseps1",
         "Volume [m$^3$]": format_number(sol_volume)},
        {"Region": "PFR",             "How it is selected": "x index $<$ ixseps1, minus core",
         "Volume [m$^3$]": format_number(pfr_volume)},
        {"Region": "First core ring", "How it is selected": "`core_edge`",
         "Volume [m$^3$]": format_number(first_core_ring_volume)},
    ]
    result.add_table(table_rows)
    return result


def check_core_source_ring_position(case):
    """
    Which rings get the core source  H(x) * H(b - x) * pressure_source,
    and does core_volume in BOUT.inp match the volume of those rings?
    """
    result = CheckResult("E. Core source ring position and core_volume",
                         "Uses the BOUT++ definition of $x$ (boutmesh.cxx, GlobalX): "
                         "it depends on the cell index only, not on $dx$.")
    result.add_formula(r"M_X = n_x - 2\,\mathrm{MXG}")
    result.add_formula(r"x_i = \frac{i + 1/2}{M_X}, \qquad i = 0, 1, \dots, M_X - 1"
                       r"\quad (i = 0 \text{ is next to the core boundary})")
    result.add_formula(r"x_{\rm sep} = \frac{\mathrm{ixseps1} - \mathrm{MXG}}{M_X}")
    result.add_formula(r"\text{ring } i \text{ is heated if } a < x_i < b"
                       r"\quad \text{for source} = H(x-a)\,H(b-x)\,p_{\rm src}")
    result.add_formula(r"\text{only core ring 1 heated:} \quad \frac{0.5}{M_X} < b < \frac{1.5}{M_X}")
    result.add_formula(r"V_{\rm heated} = \sum_{\text{heated rings}} \; \sum_{\theta \in \text{core}} \Delta V"
                       r"\;\; \overset{?}{=} \;\; \text{core\_volume}")
    result.add_formula(r"x^{\rm old}_i = \frac{\sum_{k \le i} dx_k}{\sum_k dx_k}"
                       r"\quad \text{(calculator.ipynb, not what BOUT++ uses)}")
    result.add_variable("x (size)", FROM_RESULT, "number of radial points $n_x$ (with guards)", "")
    result.add_variable("MXG", FROM_GRID, "number of x guard cells", "")
    result.add_variable("ixseps1", FROM_GRID, "x index of the first SOL ring", "")
    result.add_variable("dv", FROM_RESULT, r"cell volume $\Delta V$", "m$^3$")
    result.add_variable("dx", FROM_RESULT, "radial cell width (old method only)", "")
    for section_name in CORE_SOURCE_SECTIONS:
        result.add_variable(f"{section_name}:source", FROM_INPUT, "source expression with $H(\\cdot)$", "")
        result.add_variable(f"{section_name}:core_volume", FROM_INPUT, "volume used to spread the power", "m$^3$")
        result.add_variable(f"{section_name}:source_only_in_core", FROM_INPUT, "zero the source outside the core", "")
        result.add_variable(f"{section_name}_src", FROM_RESULT, "source in the output (cross-check)", "Pa s$^{-1}$")

    ds = case.ds
    number_of_x_guards, separatrix_index, number_of_radial_cells = get_radial_grid_sizes(ds)
    ring_table = build_radial_ring_table(ds)

    first_core_ring = ring_table.iloc[0]
    first_sol_ring = ring_table.loc[separatrix_index]
    separatrix_x = (separatrix_index - number_of_x_guards) / number_of_radial_cells

    grid_rows = [
        {"Quantity": "Radial cells without guards $M_X$", "Value": str(number_of_radial_cells)},
        {"Quantity": "MXG", "Value": str(number_of_x_guards)},
        {"Quantity": "ixseps1", "Value": str(separatrix_index)},
        {"Quantity": "First core ring: $x$ centre", "Value": format_number(first_core_ring["x centre"])},
        {"Quantity": "First core ring: $x$ faces",
         "Value": f"{format_number(first_core_ring['x inner face'])} to "
                  f"{format_number(first_core_ring['x outer face'])}"},
        {"Quantity": r"Separatrix $x_{\rm sep}$", "Value": format_number(separatrix_x)},
        {"Quantity": "First SOL ring: $x$ centre", "Value": format_number(first_sol_ring["x centre"])},
        {"Quantity": r"First SOL ring: $S_{\rm rad}$ at OMP [mm]",
         "Value": format_number(first_sol_ring["Srad OMP [mm]"])},
        {"Quantity": "Upper bound $b$ to heat only core ring 1",
         "Value": f"{format_number(0.5 / number_of_radial_cells)} < b < "
                  f"{format_number(1.5 / number_of_radial_cells)}"},
        {"Quantity": "Safest $b$ (cell face)", "Value": format_number(1.0 / number_of_radial_cells)},
    ]
    result.add_table(grid_rows, "Grid and ring positions")

    window_rows = []
    volume_rows = []
    all_source_rings = []
    window_per_section = {}
    for section_name in CORE_SOURCE_SECTIONS:
        source_expression = get_input_raw(case.input_options, f"{section_name}:source")
        if source_expression is None:
            window_rows.append({"Section": f"`[{section_name}]`", "source": "(not set)", "a": "", "b": "",
                                "Heated rings (from BOUT.inp)": "", "Rings with src > 0 (nc)": "",
                                "Status": format_status("INFO")})
            continue
        window = find_heaviside_window(source_expression)
        if window is None:
            window_rows.append({"Section": f"`[{section_name}]`", "source": f"`{source_expression}`",
                                "a": "", "b": "", "Heated rings (from BOUT.inp)": "no H(x) window",
                                "Rings with src > 0 (nc)": "", "Status": format_status("INFO")})
            continue

        lower_bound, upper_bound = window
        window_per_section[section_name] = window
        core_only = get_input_bool(case.input_options, f"{section_name}:source_only_in_core")
        predicted_rings = find_rings_in_window(ring_table, lower_bound, upper_bound, core_only)
        all_source_rings.extend(predicted_rings)

        if len(predicted_rings) == 0:
            window_status = "FAIL"      # no ring centre in the window: the source is zero
        elif len(predicted_rings) > 1:
            window_status = "WARN"      # source spread over several rings
        else:
            window_status = "PASS"

        source_variable = f"{section_name}_src"
        if source_variable in ds:
            rings_in_result = find_rings_with_positive_values(ds, source_variable)
            result_text = ring_labels_as_text(ring_table, rings_in_result)
            if rings_in_result != predicted_rings:
                window_status = "FAIL"  # result file does not agree with BOUT.inp
        else:
            result_text = "not in result file"
        result.update_status(window_status)

        window_rows.append({"Section": f"`[{section_name}]`", "source": f"`{source_expression}`",
                            "a": format_number(lower_bound), "b": format_number(upper_bound),
                            "Heated rings (from BOUT.inp)": ring_labels_as_text(ring_table, predicted_rings),
                            "Rings with src > 0 (nc)": result_text,
                            "Status": format_status(window_status)})

        heated_volume = float(ring_table.loc[predicted_rings, "ring volume"].sum())
        core_volume_input = get_input_number(case.input_options, f"{section_name}:core_volume")
        volume_rows.append(make_comparison_row(result, f"[{section_name}] $V_{{\\rm heated}}$ vs core_volume",
                                               heated_volume, core_volume_input, "m$^3$"))

    result.add_table(window_rows, "Source window per section (PASS = exactly one ring heated)")
    if len(volume_rows) > 0:
        result.add_table(volume_rows, "core_volume in BOUT.inp vs volume of the heated rings")

    ring_rows = []
    for x_index, ring in select_rings_to_show(ring_table, separatrix_index).iterrows():
        if x_index in all_source_rings:
            heated_text = "yes"
        else:
            heated_text = "no"
        ring_rows.append({"x index": str(x_index), "Ring": ring["ring"],
                          "$x$ inner face": format_number(ring["x inner face"]),
                          "$x_i$ centre": format_number(ring["x centre"]),
                          "$x$ outer face": format_number(ring["x outer face"]),
                          "old notebook $x$": format_number(ring["old notebook x"]),
                          r"$S_{\rm rad}$ OMP [mm]": format_number(ring["Srad OMP [mm]"]),
                          "Ring volume [m$^3$]": format_number(ring["ring volume"]),
                          "Heated": heated_text})
    result.add_table(ring_rows, "Rings near the core boundary and the separatrix")

    figure = plot_ring_position_figure(ds, ring_table, window_per_section, all_source_rings, separatrix_x)
    result.add_figure(figure, "Left: rings along the BOUT++ x axis and the H(x) windows. "
                              "Right: red = cells that get the core source, blue = first SOL ring.")
    return result


def check_particle_balance(case):
    """Particle inventory and rates (calculator.ipynb 'important parameters' and 'recycling vs puff')."""
    result = CheckResult("F. Particle balance (last time step)",
                         "Sums over the domain without guard cells.")
    result.add_formula(r"N_s = \sum_i n_{s,i}\,\Delta V_i, \qquad \Gamma_X = \sum_i S_{X,i}\,\Delta V_i")
    result.add_formula(r"\Gamma_{\rm target} = \frac{\Gamma_{\rm recycle}}{R}, \qquad "
                       r"\Gamma_{\rm loss} = (1 - R)\,\Gamma_{\rm target}")
    result.add_formula(r"\Gamma_{\rm net} = \Gamma_{\rm puff} + \Gamma_{\rm pump} - \Gamma_{\rm loss}"
                       r" + (\Gamma_{\rm core,n} + \Gamma_{\rm core,i})"
                       r"\;\to\; 0 \text{ at steady state}")
    result.add_formula(r"\Gamma_{\rm net,n} = \Gamma_{\rm puff} + \Gamma_{\rm pump} + \Gamma_{\rm recycle}"
                       r" + \Gamma_{\rm wall} - \Gamma_{\rm iz} - \Gamma_{\rm rec} + \Gamma_{\rm core,n}"
                       r"\;\to\; 0 \text{ at steady state (neutrals only)}")
    for variable_name, meaning, unit in [("Nd+", "ion density", "m$^{-3}$"),
                                         ("Nd", "neutral density", "m$^{-3}$"),
                                         ("Ne", "electron density", "m$^{-3}$"),
                                         ("Sd_src", "puff source", "m$^{-3}$ s$^{-1}$"),
                                         ("Sd+_iz", "ionisation source of ions (> 0)", "m$^{-3}$ s$^{-1}$"),
                                         ("Sd+_rec", "recombination source of ions (< 0)", "m$^{-3}$ s$^{-1}$"),
                                         ("Sd_target_recycle", "target recycling source", "m$^{-3}$ s$^{-1}$"),
                                         ("Sd_wall_recycle", "SOL/PFR wall recycling source", "m$^{-3}$ s$^{-1}$"),
                                         ("Sd_pump", "pump sink (negative)", "m$^{-3}$ s$^{-1}$"),
                                         (CORE_IONISING_NEUTRAL_SINK,
                                          r"core ionising neutral sink $\Gamma_{\rm core,n}$ (< 0)", "m$^{-3}$ s$^{-1}$"),
                                         (CORE_IONISING_ION_SOURCE,
                                          r"core ionising ion source $\Gamma_{\rm core,i}$ (> 0)", "m$^{-3}$ s$^{-1}$"),
                                         ("dv", r"cell volume $\Delta V$", "m$^3$")]:
        result.add_variable(variable_name, FROM_RESULT, meaning, unit)
    result.add_variable(f"{RECYCLING_ION_SECTION}:target_recycle_multiplier", FROM_INPUT,
                        "recycling coefficient $R$", "")

    ds = case.ds
    recycle_multiplier = get_input_number(case.input_options,
                                          f"{RECYCLING_ION_SECTION}:target_recycle_multiplier")

    inventory_rows = []
    for variable_name in ["Nd+", "Nd", "Ne"]:
        inventory_rows.append({"Quantity": f"$N$ of `{variable_name}`", "Variable": f"`{variable_name}`",
                               "Value": format_number(integrate_over_region(ds, variable_name)),
                               "Unit": "particles"})
    rate_names = [("Puff $\\Gamma_{\\rm puff}$", "Sd_src"),
                  ("Ionisation", "Sd+_iz"),
                  ("Recombination", "Sd+_rec"),
                  ("Target recycling $\\Gamma_{\\rm recycle}$", "Sd_target_recycle"),
                  ("Wall recycling $\\Gamma_{\\rm wall}$", "Sd_wall_recycle"),
                  ("Pump $\\Gamma_{\\rm pump}$", "Sd_pump"),
                  ("Core ionising neutral sink $\\Gamma_{\\rm core,n}$", CORE_IONISING_NEUTRAL_SINK),
                  ("Core ionising ion source $\\Gamma_{\\rm core,i}$", CORE_IONISING_ION_SOURCE)]
    rates = {}
    for label, variable_name in rate_names:
        rates[variable_name] = integrate_over_region(ds, variable_name)
        inventory_rows.append({"Quantity": label, "Variable": f"`{variable_name}`",
                               "Value": format_number(rates[variable_name]), "Unit": "s$^{-1}$"})
    result.add_table(inventory_rows, "Inventories and rates")

    # Rates missing from the result file count as zero in the balances below
    rates_or_zero = {}
    for variable_name, rate in rates.items():
        if rate is None:
            rates_or_zero[variable_name] = 0.0
        else:
            rates_or_zero[variable_name] = rate

    # Core ionising moves particles from neutrals to ions: the two terms should cancel
    core_ionising_net = (rates_or_zero[CORE_IONISING_NEUTRAL_SINK]
                         + rates_or_zero[CORE_IONISING_ION_SOURCE])

    # Neutral-only balance: the only place where core ionising does not cancel.
    # Sd+_iz > 0 removes neutrals; Sd+_rec < 0 (ion loss) adds neutrals.
    neutral_net_balance = (rates_or_zero["Sd_src"] + rates_or_zero["Sd_pump"]
                           + rates_or_zero["Sd_target_recycle"] + rates_or_zero["Sd_wall_recycle"]
                           - rates_or_zero["Sd+_iz"] - rates_or_zero["Sd+_rec"]
                           + rates_or_zero[CORE_IONISING_NEUTRAL_SINK])
    neutral_rows = [
        {"Quantity": r"Core ionising $\Gamma_{\rm core,n} + \Gamma_{\rm core,i}$ (should be 0)",
         "Value": format_number(core_ionising_net), "Unit": "s$^{-1}$"},
        {"Quantity": r"Neutral net $\Gamma_{\rm net,n}$",
         "Value": format_number(neutral_net_balance), "Unit": "s$^{-1}$"},
    ]
    if rates_or_zero["Sd+_iz"] != 0:
        neutral_rows.append({"Quantity": r"$\Gamma_{\rm net,n} / \Gamma_{\rm iz}$",
                             "Value": f"{neutral_net_balance / rates_or_zero['Sd+_iz']:.2%}", "Unit": ""})
    result.add_table(neutral_rows, "Core ionising and neutral-only balance")

    if recycle_multiplier is None or rates["Sd_target_recycle"] is None:
        result.add_text("No target recycling: target loss estimate skipped.")
        return result

    ion_flow_to_targets = rates["Sd_target_recycle"] / recycle_multiplier
    ion_loss_at_targets = (1 - recycle_multiplier) * ion_flow_to_targets
    derived_rows = [
        {"Quantity": "Recycling coefficient $R$", "Value": format_number(recycle_multiplier), "Unit": ""},
        {"Quantity": r"Ion flow to targets $\Gamma_{\rm target}$",
         "Value": format_number(ion_flow_to_targets), "Unit": "s$^{-1}$"},
        {"Quantity": r"Ions lost at targets $\Gamma_{\rm loss}$",
         "Value": format_number(ion_loss_at_targets), "Unit": "s$^{-1}$"},
    ]
    puff_rate = rates["Sd_src"]
    pump_rate = rates["Sd_pump"]
    if puff_rate is not None and puff_rate != 0:
        if pump_rate is None:
            pump_rate = 0.0
        net_balance = puff_rate + pump_rate - ion_loss_at_targets + core_ionising_net
        derived_rows.append({"Quantity": r"Net $\Gamma_{\rm net}$", "Value": format_number(net_balance),
                             "Unit": "s$^{-1}$"})
        derived_rows.append({"Quantity": r"$\Gamma_{\rm net} / \Gamma_{\rm puff}$",
                             "Value": f"{net_balance / puff_rate:.2%}", "Unit": ""})
    result.add_table(derived_rows, "Target loss and net balance")
    return result


def check_volume_integrals(case):
    """The list VOLUME_INTEGRALS_TO_SHOW from the settings."""
    result = CheckResult("G. Custom volume integrals", "Items from VOLUME_INTEGRALS_TO_SHOW.")
    result.add_formula(r"I = f \sum_{i \in \text{region}} q_i\,\Delta V_i"
                       r"\qquad (q = 1 \text{ gives the region volume})")
    result.add_variable("dv", FROM_RESULT, r"cell volume $\Delta V$", "m$^3$")
    table_rows = []
    for item in VOLUME_INTEGRALS_TO_SHOW:
        if item["variable"] is not None:
            result.add_variable(item["variable"], FROM_RESULT, f"$q$ for: {item['label']}", "")
        value = integrate_over_region(case.ds, item["variable"], item["region"], item["factor"])
        if value is None:
            value_text = "not in result file"
        else:
            value_text = format_number(value)
        table_rows.append({"Item": item["label"], "$q$": f"`{item['variable']}`",
                           "Region": f"`{item['region']}`", "$f$": format_number(item["factor"]),
                           "Value": value_text, "Unit": item["unit"]})
    result.add_table(table_rows)
    return result


def check_input_vs_result(case):
    """The list INPUT_VS_RESULT_COMPARISONS from the settings."""
    result = CheckResult("H. Custom input vs result comparisons",
                         "Items from INPUT_VS_RESULT_COMPARISONS.")
    result.add_formula(r"X_{\rm expected} = f_{\rm in} \times (\text{value in BOUT.inp})")
    result.add_formula(r"X_{\rm calculated} = f_{\rm res} \sum_{i \in \text{region}} q_i\,\Delta V_i")
    result.add_variable("dv", FROM_RESULT, r"cell volume $\Delta V$", "m$^3$")
    table_rows = []
    for item in INPUT_VS_RESULT_COMPARISONS:
        result.add_variable(item["input_key"], FROM_INPUT, item["label"], "")
        if item["result_variable"] is not None:
            result.add_variable(item["result_variable"], FROM_RESULT, item["label"], "")
        input_value = get_input_number(case.input_options, item["input_key"])
        if input_value is None:
            expected_value = None
        else:
            expected_value = item["input_factor"] * input_value
        calculated_value = integrate_over_region(case.ds, item["result_variable"],
                                                 item["result_region"], item["result_factor"])
        table_rows.append(make_comparison_row(result, item["label"], calculated_value,
                                              expected_value, item["unit"]))
    result.add_table(table_rows)
    return result


def check_template_example(case):
    """
    TEMPLATE for a new check. Copy, rename, fill in, then add it to CHECKS_TO_RUN.
    Not in CHECKS_TO_RUN by default.
    """
    result = CheckResult("X. My new check", "One line explaining what is checked")
    result.add_formula(r"N_d = \sum_i n_{d,i}\,\Delta V_i")
    result.add_variable("Nd", FROM_RESULT, "neutral density", "m$^{-3}$")
    result.add_variable("d:neutral_lmax", FROM_INPUT, "neutral mean free path limit", "m")

    neutral_lmax = get_input_number(case.input_options, "d:neutral_lmax")   # from BOUT.inp
    total_neutrals = integrate_over_region(case.ds, "Nd", "all_noguards")   # from the result file

    table_rows = []
    table_rows.append(make_comparison_row(result, "neutral_lmax", neutral_lmax, 0.1, "m"))
    table_rows.append(make_comparison_row(result, "total neutrals", total_neutrals, None, "atoms"))
    result.add_table(table_rows)
    return result


# ============================================================================
# 6. PLOTS  (each returns one CheckResult that holds figures)
# ============================================================================
def plot_one_polygon(ds, axis, variable_settings):
    """Draw one 2D polygon plot. Missing variables get a grey text instead."""
    variable_name = variable_settings["name"]
    use_logscale = variable_settings.get("logscale", False)
    use_absolute = variable_settings.get("absolute", False)

    if variable_name not in ds:
        axis.text(0.5, 0.5, f"{variable_name}\nnot in result file", ha="center",
                  va="center", color="grey", transform=axis.transAxes)
        axis.set_axis_off()
        return

    data = ds[variable_name]
    title = variable_name
    if use_absolute:
        data = np.abs(data)
        title = f"|{variable_name}|"
    if use_logscale:
        number_of_positive_cells = int((data > 0).sum())
        if number_of_positive_cells == 0:
            use_logscale = False
            title = title + " (all <= 0, linear)"
        else:
            # log of 0 is undefined: cells <= 0 get a floor value 10 times below the smallest
            # positive value, so they show the lowest colour and the grid shape stays visible
            smallest_positive_value = float(data.where(data > 0).min())
            data = data.where(data > 0, smallest_positive_value / 10)

    data.hermes.clear_guards().bout.polygon(ax=axis, cmap="Spectral_r", antialias=True, logscale=use_logscale)
    axis.set_title(title, fontsize=10)


def plot_polygon_pages(case):
    """One figure per entry in POLYGON_PAGES."""
    result = CheckResult("I. 2D polygon plots", "Settings in POLYGON_PAGES.")
    result.add_text(r"Log-scale panels: cells with value $\le 0$ get the lowest colour "
                    r"(10 times below the smallest positive value). $|q|$ = absolute value is plotted.")
    for page in POLYGON_PAGES:
        variables = page["variables"]
        for variable_settings in variables:
            result.add_variable(variable_settings["name"], FROM_RESULT, page["page_title"], "")

        number_of_columns = min(4, len(variables))
        number_of_rows = int(np.ceil(len(variables) / number_of_columns))
        figure, axes = plt.subplots(number_of_rows, number_of_columns,
                                    figsize=(4 * number_of_columns, 5.5 * number_of_rows),
                                    squeeze=False)
        for panel_number, variable_settings in enumerate(variables):
            row_index = panel_number // number_of_columns
            column_index = panel_number % number_of_columns
            plot_one_polygon(case.ds, axes[row_index, column_index], variable_settings)
        # hide unused panels
        for panel_number in range(len(variables), number_of_rows * number_of_columns):
            axes[panel_number // number_of_columns, panel_number % number_of_columns].set_axis_off()
        figure.suptitle(f"{page['page_title']}  -  {case.case_id}", fontsize=12)
        figure.tight_layout()
        result.add_figure(figure, page["page_title"])
    return result


def get_profile_data(ds, column_settings, variable_names):
    """Returns (dataframe, x_values, x_label) for one PROFILE_COLUMNS entry."""
    if column_settings["kind"] == "poloidal":
        target_first = column_settings.get("target_first", True)
        profile = get_1d_poloidal_data(ds, params=variable_names, region=column_settings["region"],
                                       sepadd=column_settings["sepadd"], target_first=target_first)
        x_values = np.abs(profile["Spar"])
        if target_first:
            x_label = "Parallel distance from target [m]"
        else:
            x_label = "Parallel distance from midplane [m]"
    else:
        profile = get_1d_radial_data(ds, params=variable_names, region=column_settings["region"])
        x_values = profile["Srad"]
        x_label = "X - X_sep [m]"
    return profile, x_values, x_label


def plot_profiles(case):
    """PROFILE_ROWS x PROFILE_COLUMNS in one figure (calculator.ipynb 'Quick check profiles')."""
    result = CheckResult("J. 1D profiles", "Settings in PROFILE_ROWS and PROFILE_COLUMNS. "
                                           "sepadd = 0 is the first SOL ring.")
    result.add_formula(r"S_\parallel = \sum_{\text{cells}} \Delta l_{\rm pol}\,\frac{|B|}{|B_{\rm pol}|}"
                       r"\quad (0 \text{ at the target})")
    result.add_formula(r"S_{{\rm rad},i} = \sum_{k \le i} dr_k - \tfrac{1}{2}\,dr_i - S_{\rm sep}"
                       r"\quad (0 \text{ at the separatrix})")
    result.add_variable("dpol, Bxy, Bpxy", FROM_RESULT, r"poloidal cell length $\Delta l_{\rm pol}$, $|B|$, $|B_{\rm pol}|$", "")
    result.add_variable("dr", FROM_RESULT, "radial cell width", "m")

    all_variable_names = []
    for row_settings in PROFILE_ROWS:
        for variable_name in row_settings["variables"]:
            if variable_name not in all_variable_names:
                all_variable_names.append(variable_name)
                result.add_variable(variable_name, FROM_RESULT, "plotted |value| in the 1D profiles", "")

    figure, axes = plt.subplots(len(PROFILE_ROWS), len(PROFILE_COLUMNS),
                                figsize=(5 * len(PROFILE_COLUMNS), 3.5 * len(PROFILE_ROWS)),
                                squeeze=False)
    for column_index, column_settings in enumerate(PROFILE_COLUMNS):
        profile, x_values, x_label = get_profile_data(case.ds, column_settings, all_variable_names)
        for row_index, row_settings in enumerate(PROFILE_ROWS):
            axis = axes[row_index, column_index]
            for variable_name in row_settings["variables"]:
                if variable_name in profile.columns:
                    axis.plot(x_values, np.abs(profile[variable_name]), label=variable_name)
            axis.set_yscale(row_settings["yscale"])
            axis.set_ylabel(row_settings["ylabel"])
            axis.set_xlabel(x_label)
            axis.grid(alpha=0.5)
            axis.legend(fontsize=8)
            if column_settings["kind"] == "radial":
                axis.axvline(0, color="black", linestyle="--", linewidth=0.8)
            if row_index == 0:
                axis.set_title(column_settings["title"])
    figure.suptitle(f"Profiles  -  {case.case_id}", fontsize=12)
    figure.tight_layout()
    result.add_figure(figure, "1D profiles")
    return result


def plot_particle_balance_history(case):
    """Particle rates versus time (calculator.ipynb 'Particle Balance' plot)."""
    result = CheckResult("K. Particle balance versus time", "Sums over the domain without guard "
                                                            "cells, for every output time.")
    result.add_formula(r"\frac{dN}{dt} = \frac{d}{dt} \sum_i (n_{e,i} + n_{d,i})\,\Delta V_i")
    result.add_formula(r"\Gamma_{\rm net}(t) = \Gamma_{\rm puff} + \Gamma_{\rm pump}"
                       r" - \frac{1-R}{R}\,\Gamma_{\rm recycle}"
                       r" + (\Gamma_{\rm core,n} + \Gamma_{\rm core,i})")
    for variable_name in ["Ne", "Nd", "Sd_src", "Sd_pump", "Sd_target_recycle",
                          CORE_IONISING_NEUTRAL_SINK, CORE_IONISING_ION_SOURCE, "dv", "t"]:
        result.add_variable(variable_name, FROM_RESULT, "see formula", "")
    result.add_variable(f"{RECYCLING_ION_SECTION}:target_recycle_multiplier", FROM_INPUT,
                        "recycling coefficient $R$", "")

    ds_all_times = case.ds_all_times
    if ds_all_times.sizes["t"] < 2:
        result.add_text("Only one time step, nothing to plot.")
        return result

    time = ds_all_times["t"]
    total_plasma = sum_over_domain_versus_time(ds_all_times, "Ne")
    total_neutral = sum_over_domain_versus_time(ds_all_times, "Nd")
    puff_rate = sum_over_domain_versus_time(ds_all_times, "Sd_src")
    pump_rate = sum_over_domain_versus_time(ds_all_times, "Sd_pump")
    target_recycle = sum_over_domain_versus_time(ds_all_times, "Sd_target_recycle")
    core_neutral_sink = sum_over_domain_versus_time(ds_all_times, CORE_IONISING_NEUTRAL_SINK)
    core_ion_source = sum_over_domain_versus_time(ds_all_times, CORE_IONISING_ION_SOURCE)
    recycle_multiplier = get_input_number(case.input_options,
                                          f"{RECYCLING_ION_SECTION}:target_recycle_multiplier")

    final_value_rows = []
    figure, axis = plt.subplots(figsize=(10, 6))
    if total_plasma is not None and total_neutral is not None:
        rate_of_change = (total_plasma + total_neutral).differentiate("t")
        axis.plot(time, rate_of_change, label="d/dt (Ne + Nd total)", lw=2, color="black")
        final_value_rows.append({"Quantity": "$dN/dt$ at last time",
                                 "Value": format_number(float(rate_of_change[-1])), "Unit": "s$^{-1}$"})

    net_balance = 0 * time
    if puff_rate is not None:
        axis.plot(time, puff_rate, label="Puff (Sd_src)", linestyle="--", lw=2)
        net_balance = net_balance + puff_rate
    if pump_rate is not None:
        axis.plot(time, pump_rate, label="Pump (Sd_pump)", linestyle="--", lw=2)
        net_balance = net_balance + pump_rate
    if target_recycle is not None:
        axis.plot(time, target_recycle, label="Target recycle", linestyle="--", lw=2)
        if recycle_multiplier is not None:
            target_loss = (1 - recycle_multiplier) / recycle_multiplier * target_recycle
            net_balance = net_balance - target_loss
    # Core ionising: neutral sink and ion source, which cancel in Ne + Nd
    if core_neutral_sink is not None:
        axis.plot(time, core_neutral_sink, label=f"Core ionising ({CORE_IONISING_NEUTRAL_SINK})",
                  linestyle="--", lw=2)
        net_balance = net_balance + core_neutral_sink
    if core_ion_source is not None:
        net_balance = net_balance + core_ion_source
    axis.plot(time, net_balance, label="Net = puff + pump - target loss + core ionising", linestyle=":", lw=3,
              color="magenta")
    final_value_rows.append({"Quantity": r"$\Gamma_{\rm net}$ at last time",
                             "Value": format_number(float(net_balance[-1])), "Unit": "s$^{-1}$"})
    result.add_table(final_value_rows)

    axis.set_xlabel("Time [s]")
    axis.set_ylabel("Particle rate [1/s]")
    axis.set_yscale("symlog", linthresh=1e18)
    axis.set_title(f"Particle balance  -  {case.case_id}")
    axis.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    axis.grid(True, alpha=0.5)
    figure.tight_layout()
    result.add_figure(figure, "Particle rates versus time")
    return result


def find_polygon_template_variable(ds):
    """Name of a variable that has the 'geometry' attribute needed by .bout.polygon()."""
    for variable_name in ["Nd", "Ne", "Td", "Te", "Pd", "Nd+"]:
        if variable_name in ds and "geometry" in ds[variable_name].attrs:
            return variable_name
    for variable_name in ds.data_vars:
        if "geometry" in ds[variable_name].attrs and set(ds[variable_name].dims) == {"x", "theta"}:
            return variable_name
    raise ValueError("No variable with a 'geometry' attribute found for the polygon plot")


def plot_ring_position_figure(ds, ring_table, window_per_section, source_rings, separatrix_x):
    """
    Left : ring cells along the BOUT++ x axis, with the H(x) window of each source.
    Right: 2D map, red = cells that get the core source, blue = first SOL ring.
    """
    number_of_x_guards, separatrix_index, number_of_radial_cells = get_radial_grid_sizes(ds)
    figure, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1.6, 1]})

    # ---- left: ring strip ----
    strip_axis = axes[0]
    largest_upper_bound = 0.0
    for window in window_per_section.values():
        if np.isfinite(window[1]):
            largest_upper_bound = max(largest_upper_bound, window[1])
    x_axis_limit = min(1.0, max(2.0 * largest_upper_bound, 8.0 / number_of_radial_cells))

    for x_index, ring in ring_table.iterrows():
        if ring["x inner face"] > x_axis_limit:
            break
        if ring["is core"]:
            ring_colour = "lightsteelblue"
        else:
            ring_colour = "navajowhite"
        strip_axis.axvspan(ring["x inner face"], ring["x outer face"], ymin=0.0, ymax=0.5,
                           facecolor=ring_colour, edgecolor="black", linewidth=0.8)
        strip_axis.plot(ring["x centre"], 0.25, "k.", markersize=6)
        strip_axis.text(ring["x centre"], 0.27, ring["ring"], ha="center", va="bottom",
                        rotation=90, fontsize=8)

    window_colours = ["red", "green", "purple", "orange"]
    for window_number, (section_name, window) in enumerate(window_per_section.items()):
        lower_bound = max(window[0], 0.0)
        upper_bound = min(window[1], x_axis_limit)
        bar_bottom = 0.55 + 0.1 * window_number
        strip_axis.axvspan(lower_bound, upper_bound, ymin=bar_bottom, ymax=bar_bottom + 0.08,
                           color=window_colours[window_number % len(window_colours)], alpha=0.6)
        strip_axis.text(upper_bound, bar_bottom + 0.04, f"  [{section_name}] H window",
                        va="center", fontsize=9)
    if separatrix_x <= x_axis_limit:
        strip_axis.axvline(separatrix_x, color="black", linestyle="--", label="separatrix")
        strip_axis.legend(loc="upper right")
    strip_axis.set_xlim(0, x_axis_limit)
    strip_axis.set_ylim(0, 1)
    strip_axis.set_yticks([])
    strip_axis.set_xlabel("BOUT++ x = (i + 0.5) / MX   (as used in BOUT.inp)")
    strip_axis.set_title("Radial rings (blue = core, orange = SOL) and source windows")

    # ---- right: 2D map ----
    # Copy a plasma/neutral variable as template: it carries the attributes ("geometry", ...)
    # that the polygon plot needs. "dv" does not have them after unnormalisation.
    template_variable = ds[find_polygon_template_variable(ds)].transpose("x", "theta")
    ring_map_values = np.zeros(template_variable.shape)
    core_theta_indices = _select_region(ds, "core_noguards")[1]
    for x_index in source_rings:
        ring_map_values[x_index, core_theta_indices] = 1.0
    domain_theta_indices = _select_region(ds, "all_noguards")[1]
    ring_map_values[separatrix_index, domain_theta_indices] = 2.0
    ring_map = template_variable.copy(data=ring_map_values)
    ring_map.name = "ring map"
    ring_map.bout.polygon(ax=axes[1], cmap=ListedColormap(["whitesmoke", "red", "blue"]),
                          vmin=0, vmax=2, antialias=True, add_colorbar=False)
    axes[1].set_title("red = core source cells, blue = first SOL ring")

    figure.tight_layout()
    return figure


# ============================================================================
# 7. WHICH CHECKS AND PLOTS TO RUN  (add your new function names here)
# ============================================================================
CHECKS_TO_RUN = [
    check_input_file_values,
    check_total_power_source,
    check_gas_puff_source,
    check_plasma_volumes,
    check_core_source_ring_position,
    check_particle_balance,
    check_volume_integrals,
    check_input_vs_result,
]

PLOTS_TO_MAKE = [
    plot_profiles,
    plot_polygon_pages,
    plot_particle_balance_history,
]


def run_check_safely(check_function, case):
    """Run one check. If it crashes, return a FAIL result with the error instead of stopping."""
    try:
        return check_function(case)
    except Exception:
        result = CheckResult(f"{check_function.__name__} (crashed)", "Python error below.")
        result.update_status("FAIL")
        error_rows = []
        for line in traceback.format_exc().splitlines():
            error_rows.append({"Python error (traceback)": "`" + line.replace("`", "'") + "`"})
        result.add_table(error_rows)
        return result


def run_checks(case, make_plots=True, make_history_plot=True):
    """Run CHECKS_TO_RUN (and PLOTS_TO_MAKE). Returns a list of CheckResult."""
    results = []
    for check_function in CHECKS_TO_RUN:
        print(f"  running {check_function.__name__} ...")
        results.append(run_check_safely(check_function, case))

    if make_plots:
        for plot_function in PLOTS_TO_MAKE:
            if plot_function is plot_particle_balance_history and not make_history_plot:
                continue
            print(f"  running {plot_function.__name__} ...")
            results.append(run_check_safely(plot_function, case))
    return results


# ============================================================================
# 8. OUTPUT  (the work is done in check_hermes_case_output.py)
# ============================================================================
def write_report(case, results, output_directory=None):
    """
    Write the Markdown + PDF report and update the version history.
    Default output_directory: REPORT_DIRECTORY/<case_id>/
    Returns (markdown_path, pdf_path). pdf_path is None if LaTeX failed.
    """
    if output_directory is None:
        output_directory = REPORT_DIRECTORY / case.case_id
    return write_report_files(case, results, output_directory, PASS_TOLERANCE, WARN_TOLERANCE)


# ============================================================================
# 9. COMMAND LINE
# ============================================================================
def main():
    global EXPECTED_PUFF_RATE, EXPECTED_PUFF_ENERGY_EV, PASS_TOLERANCE, WARN_TOLERANCE

    parser = argparse.ArgumentParser(description="Quick check of Hermes-3 cases, with a Markdown + PDF report")
    parser.add_argument("-i", "--input", nargs="+", default=[DEFAULT_CASE_ID],
                        help="Case ID(s) = folder name(s) in hms_output/")
    parser.add_argument("-o", "--output-dir", default=str(REPORT_DIRECTORY),
                        help="Folder for the reports (one sub-folder per case)")
    parser.add_argument("--no-plots", action="store_true", help="Tables only (faster)")
    parser.add_argument("--no-history", action="store_true",
                        help="Skip the particle balance vs time plot (reads all time steps)")
    parser.add_argument("--no-squash", action="store_true", help="Read BOUT.dmp.*.nc, not BOUT.squash.nc")
    parser.add_argument("--expected-puff-rate", type=float, default=EXPECTED_PUFF_RATE,
                        help="Expected puff rate [1/s], e.g. 1e21")
    parser.add_argument("--expected-puff-energy", type=float, default=EXPECTED_PUFF_ENERGY_EV,
                        help="Expected energy per puffed atom [eV], e.g. 0.02533")
    parser.add_argument("--pass-tol", type=float, default=PASS_TOLERANCE, help="PASS tolerance (0.01 = 1%%)")
    parser.add_argument("--warn-tol", type=float, default=WARN_TOLERANCE, help="WARN tolerance (0.1 = 10%%)")
    args = parser.parse_args()

    EXPECTED_PUFF_RATE = args.expected_puff_rate
    EXPECTED_PUFF_ENERGY_EV = args.expected_puff_energy
    PASS_TOLERANCE = args.pass_tol
    WARN_TOLERANCE = args.warn_tol

    plt.switch_backend("Agg")   # no windows, only files

    for case_id in args.input:
        print(f"\n##### Checking {case_id} #####")
        case = load_case(case_id, use_squash=not args.no_squash)
        results = run_checks(case, make_plots=not args.no_plots,
                             make_history_plot=not args.no_history)
        print_results(results)
        write_report(case, results, Path(args.output_dir) / case_id)


if __name__ == "__main__":
    main()

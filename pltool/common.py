# hermes3_2dplots/common.py
import xarray as xr
import matplotlib
import pandas as pd
import numpy as np
import matplotlib
import argparse
import matplotlib.pyplot as plt
import os, sys, pathlib, shlex, subprocess, glob


import xbout
import scipy
import xhermes
# from xhermes import *

sys.path.append(r"/Users/zero/workspace/phd_work/hms_output/sdtools/")
# sys.path.append(r"/users/jpm590/2dspace/post-processing/sdtools")
import hermes3

from hermes3.utils import *
from hermes3.fluxes import *
from hermes3.case_db import *
from hermes3.load import *
from hermes3.named_selections import *
from hermes3.plotting import *
from hermes3.grid_fields import *
from hermes3.accessors import *
from hermes3.selectors import *

from code_comparison.solps_pp import *
import time


# ---- Plot style config for multi-file / benchmark comparisons ----

# Hermes-3 file colours; cycles if more files than entries.
H3_COLORS = ["royalblue", "limegreen", "darkorange","magenta"]
# H3_COLORS = ["cyan", "darkorange", "magenta"]
H3_STYLES = [
    dict(ls="solid", lw=3),                
    dict(ls="solid", lw=3 ),  
    dict(ls="solid", lw=3),  
    dict(ls="solid", lw=3),  
]

# SOLPS reference curve style.
SOLPS_COLORS = ["dimgrey",  "saddlebrown", "teal", "black",]
SOLPS_STYLES = [
    dict(ls="None",        marker='.',    alpha=0.8),
    dict(ls="None",     marker='x',   alpha=0.8),
    dict(ls=(0, (1, 1)),        alpha=0.8),
    dict(ls=(0, (3, 1, 1, 1)),  alpha=0.8),
]
# keep old singular names so nothing else breaks
# SOLPS_COLOR = SOLPS_COLORS[0]
# SOLPS_STYLE = SOLPS_STYLES[0]


def build_base_parser():
    """Base parser shared by multi and single plot scripts."""
    p = argparse.ArgumentParser(
            description = "Assign a input folder and output directory"
    )

    p.add_argument(
        "-i", "--input",
        required=True,
        nargs="+",                 
        type=str,
        help="Input file IDs or patterns (e.g. 260207*)"
    )
    # If multiple files to be compared, input the legend name manually
    p.add_argument(
        "-in", "--input_name",
        required=False,
        nargs="+",                 
        type=str,
        help="Input file names for printing on plot"
    )
    p.add_argument(
        "-si", "--solps_input",
        required=False, 
        nargs="+", 
        type=str, 
        default=None,
        help="SOLPS case folder IDs or patterns (each folder must contain balance.nc)"
    )
    p.add_argument(
        "-sin", 
        "--solps_name",
        required=False, 
        nargs="+", 
        type=str, 
        default=None,
        help="Short names for SOLPS cases, used in legends"
    )
    p.add_argument("-o", "--output", 
            required=True, 
            type=str, 
            help="Path to plots folder, better to use date as note, e.g. YYMMDD"
            )
    p.add_argument("-r", "--region_rad",  
            type=str, 
            default="omp", 
            help="omp, {inner/outer}_{lower/upper}_target ... for more see doc"
            )
    p.add_argument("-p", "--region_pol",  
            type=str, 
            default="outer_lower_sol", 
            help="Must specify sepadd/sepdist ... for more see doc"
            )
    p.add_argument("--sepadd",  
            type=int, 
            default=1, 
            help="Index of the SOL ring based on nx. Default SOL ring = 1"
            )
    p.add_argument("-s", "--scale",  
            type=str, 
            default="linear", 
            help="linear or log")
    p.add_argument("-m", "--mode",
            required=True,
            nargs="+",
            type=str,
            help="benchmark? polygon? multifiles? multi_benchmark?")
    p.add_argument("--solps",
            action="store_true",
            help="In polygon mode, add a side-by-side SOLPS polygon column (Hermes-3 | SOLPS)")


    return p

def read_files(input_ids, input_name=None):

    # case_dir = "/users/jpm590/scratch/"
    case_dir = "/Users/zero/workspace/phd_work/hms_output"
    db = CaseDB(
        case_dir = case_dir,
        # grid_dir = r"/users/jpm590/gridfile"
        grid_dir = r"/Users/zero/workspace/phd_work/hms_output/gridfile/"
        # grid_dir = r"/users/jpm590/neutralimit_gcc12_9497476/hermes-3/build-mc-master"
        # grid_dir = r"/users/jpm590/benchspace/hermes-3/build-mc-master"
        # grid_dir = r"/users/jpm590/neutralrun/hermes-3/master"
        # grid_dir = r"/users/jpm590/2dspace/hermes-3/build-mc-master"
    )


    expanded_ids = []

    for pattern in input_ids:

        exact_path = os.path.join(case_dir, pattern)
        if any(c in pattern for c in ["*", "?", "["]):
            search_pattern = exact_path
            matches = glob.glob(search_pattern)
        elif os.path.exists(exact_path):
            search_pattern = exact_path
            matches = [exact_path]
        else:
            search_pattern = exact_path + "*"
            matches = glob.glob(search_pattern)

        # print("Searching:", search_pattern)
        # print("Matches:", matches)

        if matches:
            # keep only case names, not full path
            expanded_ids.extend(os.path.basename(m) for m in matches)
        else:
            print(f"No match for {pattern}, using literal")
            expanded_ids.append(pattern)

    print("Expanded IDs:", expanded_ids)


    toload = []
    for i, path in enumerate(expanded_ids):

        file_name = os.path.basename(path)

        if input_name:
            print(f"input name[{i}] = {input_name[i]}")
            toload.append(
                dict(
                    name=input_name[i],          # <-- filename here
                    id=file_name,            # or keep full path if needed
                    unnormalise_geom=True,
                    use_xhermes=True,
                    squash=True
                )
            )
        else: 
            toload.append(
                dict(
                    name=expanded_ids[i],          # <-- filename here
                    id=file_name,            # or keep full path if needed
                    unnormalise_geom=True,
                    use_xhermes=True,
                    squash=True
                )
            )

    cs = {}
    for case in toload:
        cs[case["name"]] = db.load_case_2D(case["id"], use_squash = case["squash"], verbose = True)
        cs[case["name"]].extract_2d_tokamak_geometry()
        # Poloidal projection of the parallel neutral flow, to match SOLPS fort.44
        ds = cs[case["name"]].ds
        bratio = abs(ds["Bpxy"] / abs(ds["Bxy"]))
        ds["Vd_pol"] = ds["Vd"] * bratio
        if "NVd" in ds:
            ds["NVd_pol"] = ds["NVd"] * bratio

    return cs

# def read_solps_balance():

#     # SOLPScase init includes "balance.nc", so only path is required
#     # input_path = "/Users/zero/workspace/phd_work/hms_output/SOLPS_limiter_off"
#     # input_path = "/Users/zero/workspace/phd_work/hms_output/SOLPS_Luciani_off"
#     input_path = "/Users/zero/workspace/phd_work/hms_output/SOLPS_Luciani_off_cx100/"
#     # input_path = "/Users/zero/workspace/phd_work/hms_output/SOLPS_Luciani_off_cx100_1e22/"
#     # input_path = "/users/jpm590/2dspace/post-processing/xhermes"
#     bal = SOLPScase(input_path)

#     return bal
SOLPS_DIR = "/Users/zero/workspace/phd_work/hms_output"
DEFAULT_SOLPS = ["SOLPS_Luciani_off"]

def read_solps_balance(input_ids=None, input_name=None):
    """Return {name: SOLPScase}, mirroring read_files() for Hermes-3."""
    if not input_ids:
        input_ids = DEFAULT_SOLPS

    expanded = []
    for pattern in input_ids:
        exact = os.path.join(SOLPS_DIR, pattern)
        if any(c in pattern for c in ["*", "?", "["]):
            matches = sorted(glob.glob(exact))
        elif os.path.isdir(exact):
            matches = [exact]
        else:
            matches = sorted(glob.glob(exact + "*"))
        if not matches:
            print(f"No SOLPS match for {pattern}, using literal")
            matches = [exact]
        expanded.extend(matches)

    print("Expanded SOLPS paths:", expanded)

    if input_name and len(input_name) != len(expanded):
        raise ValueError(
            f"{len(expanded)} SOLPS cases != {len(input_name)} names "
            f"(a glob pattern may have expanded to several folders)"
        )

    bals = {}
    for i, path in enumerate(expanded):
        name = input_name[i] if input_name else os.path.basename(path.rstrip("/"))
        print(f"Reading SOLPS [{i}]: {path} -> {name}")
        bals[name] = adapt_solps_conventions(SOLPScase(path))   # convention fix applied here
    return bals



def adapt_solps_conventions(bal):
    """Map solps_pp fields onto Hermes-3 benchmark names, without editing solps_pp."""
    b = bal.bal
    b["Pd"] = b["Pa"]                                             # d pressure  = atom pressure (Pa)
    b["Td"] = b["Tn"]                                             # d temperature = neutral T (Tn), not Ta
    b["Edd+_cx"] = b["eirene_mc_eapl_shi_bal"].sum(axis=2) / b["vol"]  # CX energy [W/m3]
    b["Nd"] = b["Na"]
    b["Nd+"] = b["Ne"] # quasi neutrality
    # --- neutral poloidal flow from fort.44 (atoms only, = Hermes-3 'd') ---
    try:
        b["pfluxa"] = read_fort44_field(bal, "pfluxa")
    except Exception as e:
        print(f"[warn] no pfluxa for {bal.path}: {e}")
    else:
        Na = b["Na"]
        b["Vd_pol"] = np.divide(b["pfluxa"], Na, out=np.zeros_like(Na), where=Na != 0)
        b["NVd_pol"] = Na * 2 * constants("mass_p") * b["Vd_pol"]
        # pfluxa is positive along +ix. Hermes-3's Vd_pol keeps the sign of the
        # PARALLEL Vd, so if Bpol < 0 on a ring the two are mirrored. Uncomment
        # to align (no-op where Bpol > 0):
        # b["Vd_pol"]  *= np.sign(bal.g["Bpol"])
        # b["NVd_pol"] *= np.sign(bal.g["Bpol"])
        # Bonus, free from the same file: radial neutral flow
                # Parallel = poloidal / (Bpol/Btot), i.e. assume the neutral flow that
        # carries the poloidal flux is field-aligned (Hermes-3's own assumption).
        # Signed Bpol, so the parallel sign comes out right even if Bpol < 0.
        Bpol, Btot = b["Bpol"], b["Btot"]
        bratio = np.divide(Btot, Bpol, out=np.zeros_like(Btot), where=Bpol != 0)
        b["Vd"] = b["Vd_pol"] * bratio
        b["NVd"] = b["NVd_pol"] * bratio

        rfluxa = read_fort44_field(bal, "rfluxa")
        b["Vd_perp"] = np.divide(rfluxa, Na, out=np.zeros_like(Na), where=Na != 0)
        b["NVd_perp"] = Na * 2 * constants("mass_p") * b["Vd_perp"]

    bal.params = list(b.keys())

    return bal

def read_fort44_field(bal, var_name):
    """
    Read one '*eirene data field <var>' block from <case>/fort.44 and return it
    padded with zeros to the balance grid shape (nx, ny).

    fort.44 is written without guard cells, in Fortran order (poloidal fastest).
    """
    nx, ny = bal.g["nx"], bal.g["ny"]
    path = os.path.join(bal.path, "fort.44")

    if os.path.exists(path):
        with open(path, "r") as f:
            lines = f.readlines()
        data, found = [], False
        for line in lines:
            if line.startswith(f"*eirene data field {var_name}"):
                found = True
                continue
            if found:
                if line.startswith("*"):
                    break
                data.extend(float(x) for x in line.split())
        if not found:
            raise KeyError(f"{var_name} not found in {path}")
        out = np.zeros((nx, ny))
        out[1:-1, 1:-1] = np.array(data).reshape((nx - 2, ny - 2), order="F")
        return out

    # Fallback: balance.nc carries the same fort.44 tallies (dab2-style layout,
    # 5 trailing EIRENE cells). Covers the cases that have no fort.44 on disk.
    if var_name in bal.bal:
        return bal.bal[var_name][:-5, :, 0]

    raise FileNotFoundError(f"no fort.44 and no {var_name} in balance.nc: {bal.path}")


def format_legend_and_axes(fig, ax, n_plots, xlabel):
    handles, labels = ax[0, 0].get_legend_handles_labels()
    print(handles, labels)
    
    fig.legend(handles, labels, loc='lower center', ncol=3)
    
    for i, axi in enumerate(ax.flat):
        if i >= n_plots:
            axi.axis('off')
            continue
            
        if i >= (n_plots - 3):
            axi.set_xlabel(xlabel)
            
    plt.tight_layout(rect=[0, 0.02, 1, 1])

def setup_matplotlib() -> None:
    """
    Central place for Matplotlib rcParams / styling.
    Call this once at program start (e.g. in run_single_plots / run_multi_plots).
    # Use a reasonable style without being too fancy
    try:
        plt.style.use("seaborn-v0_8-colorblind")
    except Exception:
        # Fall back silently if style not available
        pass
    """
    # Base DPI and font
    plt.rcParams["figure.dpi"] = 120
    plt.rcParams["savefig.dpi"] = 300
    plt.rcParams["font.size"] = 14

    # Axes and lines
    plt.rcParams["axes.grid"] = False  # you turn grid on per-axis
    plt.rcParams["axes.labelsize"] = 14
    plt.rcParams["axes.titlesize"] = 14
    plt.rcParams["lines.linewidth"] = 1.5
    # If compare multiple files
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=H3_COLORS)
    # plt.rcParams["axes.prop_cycle"] = plt.cycler(color=plt.cm.gist_rainbow.colors)

    # Legend
    plt.rcParams["legend.frameon"] = True  # do you want the frame or not
    plt.rcParams["legend.fontsize"] =  14


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
from xhermes import *

sys.path.append(r"/users/jpm590/2dspace/post-processing/sdtools")
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


def build_base_parser():
# def build_base_parser(description: str) -> argparse.ArgumentParser:
    """Base parser shared by multi and single plot scripts."""
    p = argparse.ArgumentParser(
            description = "Assign a input folder and output directory"
    )
    # p.add_argument("-i", "--input", required=True, type=str, help="Name of netcdf files folder, not path")
    p.add_argument(
        "-i", "--input",
        required=True,
        nargs="+",                 # <-- key change
        type=str,
        help="Input file IDs or patterns (e.g. 260207*)"
    )
    p.add_argument("-o", "--output", required=True, type=str, help="Path to plots folder, better to use date as note, e.g. YYMMDD")
    p.add_argument("-r", "--region_rad",  type=str, default="omp", help="omp, {inner/outer}_{lower/upper}_target ... for more see doc")
    p.add_argument("-p", "--region_pol",  type=str, default="outer_lower", help="Must specify sepadd/sepdist ... for more see doc")
    p.add_argument("--sepadd",  type=int, default=1, help="Index of the SOL ring based on nx. Default SOL ring = 1")
    p.add_argument("-s", "--scale",  type=str, default="linear", help="linear or log")
    # p.add_argument("--isolps", type=str, help="Path to SOLPS's balance.nc")

    return p

def read_files(input_ids):

    case_dir = "/users/jpm590/scratch/"
    db = CaseDB(
        case_dir = case_dir,
        grid_dir = r"/users/jpm590/neutralimit_gcc12_9497476/hermes-3/build-mc-master"
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

        print("Searching:", search_pattern)
        print("Matches:", matches)

        if matches:
            # keep only case names, not full path
            expanded_ids.extend(os.path.basename(m) for m in matches)
        else:
            print(f"No match for {pattern}, using literal")
            expanded_ids.append(pattern)

    print("Expanded IDs:", expanded_ids)


    # Reduce common prefix, get the short names for plots
    prefix = os.path.commonprefix(expanded_ids)
    short_names = [name[len(prefix):] for name in expanded_ids]
    print(short_names)

    # Build toload
    toload = []
    for i, path in enumerate(expanded_ids):

        file_name = os.path.basename(path)

        toload.append(
            dict(
                name=short_names[i],          # <-- filename here
                id=file_name,            # or keep full path if needed
                unnormalise_geom=True,
                use_xhermes=True,
                squash=True
            )
        )

        print(file_name)

    cs = {}
    for case in toload:
        cs[case["name"]] = db.load_case_2D(case["id"], use_squash = case["squash"], verbose = True)
        cs[case["name"]].extract_2d_tokamak_geometry()

    return cs

def read_solps_balance():

    # SOLPScase init includes "balance.nc", so only path is required
    input_path = "/users/jpm590/2dspace/post-processing/xhermes"
    bal = SOLPScase(input_path)

    return bal


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
    plt.rcParams["font.size"] = 13

    # Axes and lines
    plt.rcParams["axes.grid"] = False  # you turn grid on per-axis
    plt.rcParams["axes.labelsize"] = 13
    plt.rcParams["axes.titlesize"] = 13
    plt.rcParams["lines.linewidth"] = 1.5

    # Legend
    plt.rcParams["legend.frameon"] = True  # do you want the frame or not
    plt.rcParams["legend.fontsize"] =  9


# hermes3_2dplots/common.py
import xarray as xr
import matplotlib
import pandas as pd
import numpy as np
import matplotlib
import argparse
import matplotlib.pyplot as plt
import os, sys, pathlib, shlex, subprocess

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

import time


def build_base_parser():
# def build_base_parser(description: str) -> argparse.ArgumentParser:
    """Base parser shared by multi and single plot scripts."""
    p = argparse.ArgumentParser(
            description = "Assign a input folder and output directory"
    )
    p.add_argument("-i", "--input", required=True, type=str, help="Name of netcdf files folder, not path")
    p.add_argument("-o", "--output", required=True, type=str, help="Path to plots folder, better to use date as note, e.g. YYMMDD")
    p.add_argument("-r", "--region_rad",  type=str, default="omp", help="omp, {inner/outer}_{lower/upper}_target ... for more see doc")
    p.add_argument("-p", "--region_pol",  type=str, default="outer_lower", help="Must specify sepadd/sepdist ... for more see doc")
    p.add_argument("--sepadd",  type=int, default=1, help="Index of the SOL ring based on nx. Default SOL ring = 1")
    p.add_argument("-s", "--scale",  type=str, default="linear", help="linear or log")

    return p

def read_files(input_id):

    db = CaseDB(
        case_dir = r"/users/jpm590/scratch/",
        grid_dir = r"/users/jpm590/2dspace/hermes-3/build-mc-master"
    )
    
    toload = [
        dict(name="MAST-U", id=input_id , unnormalise_geom = True, use_xhermes = True, squash = True)

    ]
    cs = {}
    for case in toload:
        cs[case["name"]] = db.load_case_2D(case["id"], use_squash = case["squash"], verbose = True)
        cs[case["name"]].extract_2d_tokamak_geometry()


    return cs


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
    plt.rcParams["font.size"] = 11

    # Axes and lines
    plt.rcParams["axes.grid"] = False  # you turn grid on per-axis
    plt.rcParams["axes.labelsize"] = 11
    plt.rcParams["axes.titlesize"] = 11
    plt.rcParams["lines.linewidth"] = 1.5

    # Legend
    plt.rcParams["legend.frameon"] = True  # do you want the frame or not
    plt.rcParams["legend.fontsize"] =  9


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


p = argparse.ArgumentParser(
        description = "Assign a input folder and output directory"
)
p.add_argument("-i", "--input", required=True, type=str, help="Name of netcdf files folder, not path")
p.add_argument("-o", "--output", required=True, type=str, help="Path to plots folder, better to use date as note, e.g. YYMMDD")
p.add_argument("-r", "--region_rad",  type=str, default="omp", help="omp, {inner/outer}_{lower/upper}_target ... for more see doc")
p.add_argument("-p", "--region_pol",  type=str, default="outer_lower", help="Must specify sepadd/sepdist ... for more see doc")
p.add_argument("--sepadd",  type=int, default=1, help="Index of the SOL ring based on nx. Default SOL ring = 1")
p.add_argument("-s", "--scale",  type=str, default="linear", help="linear or log")


args = p.parse_args()

def read_file(input_id):

    # Read file 

    db = CaseDB(
        case_dir = r"/users/jpm590/scratch/",
        #case_dir = r"/users/jpm590/2dspace/run/",
        grid_dir = r"/users/jpm590/2dspace/hermes-3/build-mc-master"
    )
    
    toload = [
        # dict(name="MAST-U", id="250929-shorter-test", unnormalise_geom = True, use_xhermes = True, squash = True)
        # dict(name="MAST-U", id="251007-2D-MASTU", unnormalise_geom = True, use_xhermes = True, squash = True)
        dict(name="MAST-U", id=input_id , unnormalise_geom = True, use_xhermes = True, squash = True)
        # dict(name="MAST-U", id="251107-tuned-puff-1e21", unnormalise_geom = True, use_xhermes = True, squash = True)

    ]
    cs = {}
    for case in toload:
        cs[case["name"]] = db.load_case_2D(case["id"], use_squash = case["squash"], verbose = True)
        cs[case["name"]].extract_2d_tokamak_geometry()

    print("Loaded file")
    return cs


def plot_multi_profiles_fieldline(cs, region_pol, idx_ring_array):
    plots = [
        ("Te",      "e temperature",        "T [eV]",                    False),
        ("Td+",     "d+ temperature",       "T [eV]",                    False),
        ("Td",      "d temperature",        "T [eV]",                    False),
        ("Ne",      "e density",            "density [m^-3]",            True),
        ("Nd+",     "d+ density",           "density [m^-3]",            True),
        ("Nd",      "d density",            "density [m^-3]",            True),
        ("Sd+_rec", "Recombination",        "Rate [m^-3 s^-1]",          True),
        ("Sd+_iz",  "Ionisation",           "Rate [m^-3 s^-1]",          True),
        ("Edd+_cx", "Charge exchange",      "Energy [W m^3]",   True),
        ("Fdd+_cx", "Charge exchange",      "Momentum [kg m^-2 s^-2]",   True),
        ("NVd",     "d Parallel momentum",  "Momentum [kg m^-2 s^-2]",   True),
        ("NVd+",    "d+ Parallel momentum", "Momentum [kg m^-2 s^-2]",   True),
    ]
   
    ncols = len(plots)
    
    fig, ax = plt.subplots(int(ncols/3), 3, figsize= (13.5, 18), squeeze=False)

    ds = cs["MAST-U"].ds.isel(t=-1)
    xpt_spar_list = []
    for ring in idx_ring_array:
        print(f"Plotting {ring} ....")
        df = get_1d_poloidal_data(ds, params=[p[0] for p in plots] + ["Bpxy"], region=region_pol, sepadd=ring)
    
        for idx, (param, title, ylabel, logy) in enumerate(plots):
            r, c = divmod(idx, 3)
            axi = ax[r,c]
       
            # B_values = np.array(df["Bpxy"][9:15])
            # min_B = np.min(B_values)
            # min_B_loc = np.argmin(B_values)
            # print(f"X-point is at {min_B_loc} with value {min_B}")
            N = len(df["Bpxy"].values)
            search_start = int(N * 0.2)
            search_end   = int(N * 0.8)
            
            B_mid = df["Bpxy"][search_start:search_end]
            
            local_min = np.argmin(B_mid)
            xpt_index = search_start + local_min
            xpt_value = df["Bpxy"][xpt_index]
            xpt_spar  = df["Spar"][xpt_index]
            xpt_spar_list.append(float(xpt_spar))
            
            print(f"X-point index = {xpt_index}")
            print(f"X-point Bpxy = {xpt_value}")
            print(f"X-point Spar = {xpt_spar}")
        
            # print(np.min(B_value))
            axi.plot(np.abs(df["Spar"][::-1]), np.abs(df[param]), label=f"ring = {ring}")
            axi.set_title(title)
            axi.set_xlabel("Spar [m]")
            axi.set_ylabel(ylabel)
            # axi.axvline(xpt_spar, color='r', linestyle='--', alpha = 0.5)
            axi.grid(True, alpha=0.5)
        
            if logy:
                axi.set_yscale("log")

        xpt_min = min(xpt_spar_list)
        xpt_max = max(xpt_spar_list)
        print(f"X-point Spar band: {xpt_min} → {xpt_max}")
        
        for axi in ax.flat:
            axi.legend()
            axi.axvspan(xpt_min, xpt_max, color='red', alpha=0.2)

    plt.tight_layout()
    fig.savefig("multi_rings_profiles.png")

def main():

    print(args.scale)
    case = read_file(args.input)

    '''For test'''
    # make_plot_spar(case)

    '''
    Plot all function

    Usage:
        region_rad = "omp"
        region_pol = "outer_lower"
        idx_ring = 15 # deps on nx, e.g. 0 - 19 

    '''
    sepadd_array = [0, 1, 2, 3, 4]
    # plot_single_profiles(case, args.region_rad, args.region_pol, args.sepadd)
    plot_multi_profiles_fieldline(case, args.region_pol, sepadd_array)
if __name__ == "__main__":
    
    start_time = time.time()
    # if not os.path.exists(f"./{figures_pdf_path}"):
    #     os.makedirs(figures_pdf_path)
    # if not os.path.exists(f"./{figures_png_path}"):
    #     os.makedirs(figures_png_path)
    main()
    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds")
    

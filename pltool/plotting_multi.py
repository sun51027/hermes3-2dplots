from .common import build_base_parser, read_files, setup_matplotlib
import numpy as np
import matplotlib
import argparse
import matplotlib.pyplot as plt
import os, sys
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

def plot_multi_profiles_fieldline(cs, region_pol, idx_ring_array, figures_png_path):
    plots = [
        ("Te",      "e temperature",        "T [eV]",                    False),
        ("Td+",     "d+ temperature",       "T [eV]",                    False),
        ("Td",      "d temperature",        "T [eV]",                    False),
        ("Ne",      "e density",            "density [$m^{-3}$]",            True),
        ("Nd+",     "d+ density",           "density [$m^{-3}$]",            True),
        ("Nd",      "d density",            "density [$m^{-3}$]",            True),
        ("Sd+_rec", "Recombination",        "Rate [$m^{-3} s^{-1}$]",          True),
        ("Sd+_iz",  "Ionisation",           "Rate [$m^{-3} s^{-1}$]",          True),
        ("Edd+_cx", "Charge exchange",      "Energy [$W m^{3}$]",   True),
        ("Fdd+_cx", "Charge exchange",      "Momentum [$kg m^{-2} s^{-2}$]",   True),
        ("NVd",     "d Parallel momentum",  "Momentum [$kg m^{-2} s^{-2}$]",   True),
        ("NVd+",    "d+ Parallel momentum", "Momentum [$kg m^{-2} s^{-2}$]",   True),
    ]
   
    ncols = len(plots)
    print(f"ncols = {ncols}")
    
    fig, ax = plt.subplots(int(ncols/3), 3, figsize= (10, 10), squeeze=False)

    ds = cs["MAST-U"].ds.isel(t=-1)
    xpt_spar_list = []
    for ring in idx_ring_array:
        print(f"Plotting {ring} ....")
        df = get_1d_poloidal_data(ds, params=[p[0] for p in plots] + ["Bpxy"], region=region_pol, sepadd=ring)


        ### Find X-point 

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
    
        for idx, (param, title, ylabel, logy) in enumerate(plots):
            r, c = divmod(idx, 3)
            axi = ax[r,c]
            
            axi.plot(np.abs(df["Spar"][::-1]), np.abs(df[param]), label=f"ring = {ring}")
            axi.set_title(title)
            axi.set_ylabel(ylabel)
            axi.set_xlabel("")
            axi.grid(True, alpha=0.7)
        
            if logy:
                axi.set_yscale("log")


    ## Print out the X-point band
    xpt_min = min(xpt_spar_list)
    xpt_max = max(xpt_spar_list)
    print(f"X-point Spar band: {xpt_min} → {xpt_max}")
    
    for i, axi in enumerate(ax.flat):
        axi.legend()
        axi.axvspan(xpt_min, xpt_max, color='red', alpha=0.2)
        
        if i == ncols-1 or i == ncols-2 or i == ncols-3:
            axi.set_xlabel("$S_{\\parallel}$")

    plt.tight_layout()
    fig.savefig(f"{figures_png_path}/multi_rings_profiles.png")

def run_multi_plots():

    setup_matplotlib()
    parser = build_base_parser()
    args = parser.parse_args()
    case = read_files(args.input)

    ## create output directory
    figures_png_path = args.output + "_figures_png"
    if not os.path.exists(f"./{figures_png_path}"):
        os.makedirs(figures_png_path)

    sepadd_array = [0, 1, 2, 3, 4]
    plot_multi_profiles_fieldline(case, args.region_pol, sepadd_array, figures_png_path)

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

def plot_multifiles_profiles_fieldline(cs, region_pol, ring, figures_png_path):

    # Plot sol ring = 1 ~ 5 on the same plot for each species separately
    # With x-point for each ring

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
    for case_name, case_obj in cs.items():
        print(case_name)

        ds = case_obj.ds.isel(t=-1)
        df = get_1d_poloidal_data(ds, params=[p[0] for p in plots] + ["Bpxy"], region=region_pol, sepadd=ring, target_first=True)
       

        for idx, (param, title, ylabel, logy) in enumerate(plots):
            r, c = divmod(idx, 3)
            axi = ax[r,c]
            
            axi.plot(np.abs(df["Spar"]), np.abs(df[param]), label=f"{case_name}")
            axi.set_title(title)
            axi.set_ylabel(ylabel)
            axi.set_xlabel("")
            axi.grid(True, alpha=0.7)
            if logy:
                axi.set_yscale("log")
    

    
    for i, axi in enumerate(ax.flat):
        axi.legend()
        
        if i == ncols-1 or i == ncols-2 or i == ncols-3:
            axi.set_xlabel("$S_{\\parallel}$")

    plt.tight_layout()
    fig.savefig(f"{figures_png_path}/multifiles_rings_profiles.png")

def plot_multifiles_profiles_radial(cs, region_rad, figures_png_path):
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
    
    fig, ax = plt.subplots(int(ncols/3), 3, figsize= (10, 10), squeeze=False)

    
    for case_name, case_obj in cs.items():
        print(case_name)
        ds = case_obj.ds.isel(t=-1)

        for idx, (param, title, ylabel, logy) in enumerate(plots):

            df = get_1d_radial_data(ds, params=[p[0] for p in plots], region=region_rad)
            r, c = divmod(idx, 3)
            axi = ax[r,c]
            
            axi.plot(df["Srad"], np.abs(df[param]), label=f"{case_name}")
            axi.set_title(title)
            axi.set_ylabel(ylabel)
            axi.set_xlabel("")
            axi.grid(True, alpha=0.7)
        
            if logy:
                axi.set_yscale("log")


    
    for i, axi in enumerate(ax.flat):
        axi.legend()
        
        if i == ncols-1 or i == ncols-2 or i == ncols-3:
            axi.set_xlabel("$X-X_{sep}$ [m]")

    plt.tight_layout()
    fig.savefig(f"{figures_png_path}/multifiles_radial_profiles.png")
def run_multifiles_plots():

    setup_matplotlib()
    parser = build_base_parser()
    args = parser.parse_args()
    case= read_files(args.input)

    ## create output directory
    figures_png_path = args.output + "_figures_png"
    if not os.path.exists(f"./{figures_png_path}"):
        os.makedirs(figures_png_path)

    # sepadd_array = [0, 1, 2, 3, 4]
    sepadd_array = 1
    plot_multifiles_profiles_fieldline(case, args.region_pol, sepadd_array, figures_png_path)
    plot_multifiles_profiles_radial(case, args.region_rad, figures_png_path)

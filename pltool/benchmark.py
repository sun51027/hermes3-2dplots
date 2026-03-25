from .common import build_base_parser, read_files, setup_matplotlib, read_solps_balance

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

def plot_bm_profiles_fieldline(cs, bal, region_pol, idx_ring_array, figures_png_path):

    # Plot sol ring = 1 ~ 5 on the same plot for each species separately
    # With x-point for each ring

    plots = [
        ("Te",      "e temperature",        "T [eV]",                    False),
        ("Td+",     "d+ temperature",       "T [eV]",                    False),
        ("Td",      "d temperature",        "T [eV]",                    False),
        ("Ne",      "e density",            "density [$m^{-3}$]",            True),
        ("Nd+",     "d+ density",           "density [$m^{-3}$]",            True),
        ("Nd",      "d density",            "density [$m^{-3}$]",            True),
        ("Pe",      "e pressure",            "pressure [Pa]",            True),
        ("Pd+",     "d+ pressure",           "pressure [Pa]",            True),
        ("Pd",      "d pressure",            "pressure [Pa]",            True),
        ]
    ncols = len(plots)
    fig, ax = plt.subplots(int(ncols/3), 3, figsize= (10, 10), squeeze=False)

    case_name = list(cs.keys())[0]
    ds = cs[case_name].ds.isel(t=-1)
    xpt_spar_hlist = []
    xpt_spar_slist = []
    for ring in idx_ring_array:
        print(f"Plotting {ring} ....")
        # df = get_1d_poloidal_data(ds, params=[p[0] for p in plots], region=region_pol, sepadd=ring)
        # df_solps = bal.get_1d_poloidal_data(params=[p[0] for p in plots], region=region_pol, sepadd=ring)
        df = get_1d_poloidal_data(ds, params=[p[0] for p in plots], region=region_pol, sepadd=ring, target_first=True)
        df_solps = bal.get_1d_poloidal_data(params=[p[0] for p in plots], region=region_pol, sepadd=ring, target_first=True)

        local_Rmin_h = np.argmin(df["R"].values) # return location
        print(f"Hermes location of Rmin: {local_Rmin_h}, value: {np.min(df['R'].values)}")
        local_Rmin_solps = np.argmin(df_solps["R"].values) # return location
        print(f"SOLPS  location of Rmin: {local_Rmin_solps}, value: {np.min(df_solps['R'].values)}")

        xpt_spar_hlist.append(float(df["Spar"][local_Rmin_h]))
        xpt_spar_slist.append(float(df_solps["Spar"][local_Rmin_solps]))
        print(f"The coordinate is: {df["Spar"][local_Rmin_h]}")

        for idx, (param, title, ylabel, logy) in enumerate(plots):
            r, c = divmod(idx, 3)
            axi = ax[r,c]

            # Hermes-3 
            # axi.plot(np.abs(df["Spar"][::-1]), np.abs(df[param]), label=f"Hermes-3")
            axi.plot(np.abs(df["Spar"]), np.abs(df[param]), label=f"Hermes-3")

            # SOLPS
            axi.plot(df_solps["Spar"], np.abs(df_solps[param]), label=f"SOLPS", ls = "--")
            # axi.plot(df_solps["Spar"][::-1], np.abs(df_solps[param]), label=f"SOLPS", ls = "--")

            axi.set_title(title)
            axi.set_ylabel(ylabel)
            axi.set_xlabel("")
            axi.grid(True, alpha=0.7)
        
            if logy:
                axi.set_yscale("log")

    xpt_h = min(xpt_spar_hlist)
    xpt_s = min(xpt_spar_hlist)
    print(xpt_h)
    
    for i, axi in enumerate(ax.flat):
        axi.legend()
        axi.axvline(xpt_h, color='r', alpha = 0.2)
        axi.axvline(xpt_s, color='r', alpha = 0.2, ls = '--')
        
        if i == ncols-1 or i == ncols-2 or i == ncols-3:
            axi.set_xlabel("$S_{\\parallel}$")

    plt.tight_layout()
    fig.savefig(f"{figures_png_path}/benchmark_rings.png")



def plot_bm_profiles_radial(cs, bal, region_rad, figures_png_path):
    plots = [
        ("Te",      "e temperature",        "T [eV]",                    False),
        ("Td+",     "d+ temperature",       "T [eV]",                    False),
        ("Td",      "d temperature",        "T [eV]",                    False),
        ("Ne",      "e density",            "density [$m^{-3}$]",            True),
        ("Nd+",     "d+ density",           "density [$m^{-3}$]",            True),
        ("Nd",      "d density",            "density [$m^{-3}$]",            True),
        ("Pe",      "e pressure",            "pressure [Pa]",            True),
        ("Pd+",     "d+ pressure",           "pressure [Pa]",            True),
        ("Pd",      "d pressure",            "pressure [Pa]",            True),
        # ("Sd+_rec", "Recombination",        "Rate [$m^{-3} s^{-1}$]",          True),
        # ("Sd+_iz",  "Ionisation",           "Rate [$m^{-3} s^{-1}$]",          True),
        # ("Edd+_cx", "Charge exchange",      "Energy [$W m^{3}$]",   True),
        # ("Fdd+_cx", "Charge exchange",      "Momentum [$kg m^{-2} s^{-2}$]",   True),
        # ("NVd",     "d Parallel momentum",  "Momentum [$kg m^{-2} s^{-2}$]",   True),
        # ("NVd+",    "d+ Parallel momentum", "Momentum [$kg m^{-2} s^{-2}$]",   True),
    ]
   
    ncols = len(plots)
    print(f"ncols = {ncols}")
    
    fig, ax = plt.subplots(int(ncols/3), 3, figsize= (10, 10), squeeze=False)

    case_name = list(cs.keys())[0]
    ds = cs[case_name].ds.isel(t=-1)
    
    for idx, (param, title, ylabel, logy) in enumerate(plots):
        df = get_1d_radial_data(ds, params=[p[0] for p in plots], region=region_rad)
        df_solps = bal.get_1d_radial_data(params=[p[0] for p in plots], region=region_rad)

        r, c = divmod(idx, 3)
        axi = ax[r,c]
        # Hermes-3 
        axi.plot(df["Srad"], np.abs(df[param]), label=f"Hermes-3")
        # SOLPS dist = Srad
        axi.plot(df_solps["dist"], np.abs(df_solps[param]), label=f"SOLPS", ls = "--")

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
    fig.savefig(f"{figures_png_path}/benchmark_radial.png")

def plot_bm_plasma_overlap(cs, bal, region_pol, region_rad, figures_png_path):

    #########
    #
    #   ring 1  radial
    #   T all   T all
    #   N all   N all
    #   P all   P all
    #
    ###########

    fig, ax = plt.subplots(5, 2, figsize=(10, 16), squeeze=False)

    case_name = list(cs.keys())[0]
    ds = cs[case_name].ds.isel(t=-1)

    params = ["Td", "Te", "Td+", "Nd", "Ne", "Nd+", "Pd", "Pe", "Pd+","Sd+_rec", "Sd+_iz"]
    
    # Get Poloidal Data
    df_herm_pol = get_1d_poloidal_data(ds, params=params, region=region_pol, sepadd=1, target_first=True)
    df_solps_pol = bal.get_1d_poloidal_data(params=params, region=region_pol, sepadd=1,target_first=True)

    # Get Radial Data
    df_herm_rad = get_1d_radial_data(ds, params=params, region=region_rad)
    df_solps_rad = bal.get_1d_radial_data(params=params, region=region_rad)
    
    plot_groups = [
        {"vars": ["Te", "Td+", "Td"], "ylabel": "Temperature (eV)", "yscale": "linear", "title_prefix": ""},
        {"vars": ["Ne", "Nd+", "Nd"], "ylabel": "Density (m^-3)", "yscale": "log", "title_prefix": ""},
        {"vars": ["Pe", "Pd+", "Pd"], "ylabel": "Pressure (Pa)", "yscale": "log", "title_prefix": ""},
        {"vars": ["Sd+_rec"], "ylabel": "Recombination (m^-3 s^-1)", "yscale": "log", "title_prefix": ""},
        {"vars": ["Sd+_iz"], "ylabel": "Ionisation (m^-3 s^-1)", "yscale": "log", "title_prefix": ""},

    ]

    for row, group in enumerate(plot_groups):
        for var in group["vars"]:
            # --- Poloidal plot (Col 0) ---
            # p_pol, = ax[row, 0].plot(np.abs(df_herm_pol["Spar"][::-1]), np.abs(df_herm_pol[var]), label=f"{var} (Hermes)")
            # color = p_pol.get_color()
            # ax[row, 0].plot(np.abs(df_solps_pol["Spar"][::-1]), np.abs(df_solps_pol[var]), label=f"{var} (SOLPS)", color=color, ls=":")
            p_pol, = ax[row, 0].plot(np.abs(df_herm_pol["Spar"]), np.abs(df_herm_pol[var]), label=f"{var} (Hermes)")
            color = p_pol.get_color()
            ax[row, 0].plot(np.abs(df_solps_pol["Spar"]), np.abs(df_solps_pol[var]), label=f"{var} (SOLPS)", color=color, ls=":")
            
            # --- Radial plot (Col 1) ---
            p_rad, = ax[row, 1].plot(df_herm_rad["Srad"], np.abs(df_herm_rad[var]), label=f"{var} (Hermes)")
            color_rad = p_rad.get_color()
            ax[row, 1].plot(df_solps_rad["dist"], np.abs(df_solps_rad[var]), label=f"{var} (SOLPS)", color=color_rad, ls=":")

        # Configure Poloidal axis
        ax[row, 0].set_ylabel(group["ylabel"])
        ax[row, 0].set_yscale(group["yscale"])
        if row == 2:
            ax[row, 0].set_xlabel("Parallel distance [m]")
        if row == 0:
            ax[row, 0].set_title(group["title_prefix"])
            
        # Configure Radial axis
        ax[row, 1].set_ylabel(group["ylabel"])
        ax[row, 1].set_yscale(group["yscale"])
        if row == 2:
            ax[row, 1].set_xlabel("$X-X_{sep}$ [m]")
        if row == 0:
            ax[row, 1].set_title(group["title_prefix"])
            
    # Dummy lines for linestyle legends
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='k', ls='-', label='Hermes-3'),
        Line2D([0], [0], color='k', ls='--', label='SOLPS')
    ]

    for axi in ax.flat:
        # Create legend with variables only and the line styles
        handles, labels = axi.get_legend_handles_labels()
        # Keep only the first few handles for variables (exclude duplicates from solps if we want, but here they have different labels)
        axi.legend(fontsize='small')
        axi.grid(True, alpha=0.7)

    plt.tight_layout()
    fig.savefig(f"{figures_png_path}/benchmark_all_overlap.png")

def run_benchmark():

    setup_matplotlib()
    parser = build_base_parser()
    args = parser.parse_args()
    case = read_files(args.input)
    balance = read_solps_balance()

    ## create output directory
    figures_png_path = args.output + "_figures_png"
    if not os.path.exists(f"./{figures_png_path}"):
        os.makedirs(figures_png_path)

    sepadd_array = [1]
    # plot_bm_profiles_radial(case, balance, args.region_rad, figures_png_path)
    plot_bm_profiles_fieldline(case, balance, args.region_pol, sepadd_array, figures_png_path)
    # plot_bm_plasma_overlap(case, balance, args.region_pol, args.region_rad, figures_png_path)

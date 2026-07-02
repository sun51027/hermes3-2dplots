from .common import build_base_parser, read_files, setup_matplotlib, read_solps_balance, format_legend_and_axes
import math
import numpy as np
import matplotlib
import argparse
import matplotlib.pyplot as plt
import os, sys
import xbout
import scipy
import xhermes
# from xhermes import *

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

def adapt_solps_conventions(bal):
    """Map solps_pp fields onto Hermes-3 benchmark names, without editing solps_pp."""
    b = bal.bal
    b["Pd"] = b["Pa"]                                             # d pressure  = atom pressure (Pa)
    b["Td"] = b["Tn"]                                             # d temperature = neutral T (Tn), not Ta
    b["Edd+_cx"] = b["eirene_mc_eapl_shi_bal"].sum(axis=2) / b["vol"]  # CX energy [W/m3]
    b["Nd"] = b["Na"]
    b["Nd+"] = b["Ne"] # quasi neutrality
    bal.params = list(b.keys())
    return bal


def plot_bm_profiles_fieldline(cs, bal, region_pol, idx_ring_array, figures_png_path):


    plots = {
        "state_var": [
            ("Te",  "e temperature", "T [eV]",            False),
            ("Td+", "d+ temperature","T [eV]",            False),
            ("Td",  "d temperature", "T [eV]",            False),
            ("Ne",  "e density",     "density [$m^{-3}$]", True),
            ("Nd+", "d+ density",    "density [$m^{-3}$]", True),
            ("Nd",  "d density",     "density [$m^{-3}$]", True),
            ("Pe",  "e pressure",    "pressure [Pa]",      True),
            ("Pd+", "d+ pressure",   "pressure [Pa]",      True),
            ("Pd",  "d pressure",    "pressure [Pa]",      True),
        ],
        "reactions": [
            ("Sd+_rec", "Recombination", "Rate [$m^{-3} s^{-1}$]", True),
            ("Sd+_iz",  "Ionisation",           "Rate [$m^{-3} s^{-1}$]",          True),
            ("Edd+_cx", "Charge exchange",      "Energy [$W m^{3}$]",   True),
            # ("Fdd+_cx", "Charge exchange",      "Momentum [$kg m^{-2} s^{-2}$]",   True),
            # ("NVd",     "d Parallel momentum",  "Momentum [$kg m^{-2} s^{-2}$]",   True),
            ("NVd+",    "d+ Parallel momentum", "Momentum [$kg m^{-2} s^{-2}$]",   True),
        ]
    }

    # case_name = list(cs.keys())[0]
    # ds = cs[case_name].ds.isel(t=-1)
    xpt_spar_hlist = []
    xpt_spar_slist = []

    for ring in idx_ring_array:

        for category, items in plots.items():
            nitems = len(items)
            print(f"Number of plots = {nitems} in category {category}")
            cols = 3
            rows = math.ceil(nitems / cols)
            fig, ax = plt.subplots(rows, cols, figsize=(12, 4 * rows), squeeze=False)

            for idx, (param, title, ylabel, logy) in enumerate(items):
                r, c = divmod(idx, 3)
                axi = ax[r,c]
    
                # SOLPS
                df_solps = bal.get_1d_poloidal_data(params=[p[0] for p in items], region=region_pol, sepadd=ring, target_first=True, interpolate_midplane=False, interpolate_radial=False )
                axi.plot(df_solps["Spar"], np.abs(df_solps[param]), label=f"SOLPS", ls = "--")
    
                # Find SOLPS's x-point
                local_Rmin_solps = np.argmin(df_solps["R"].values) # return location
                # print(f"SOLPS  location of Rmin: {local_Rmin_solps}, value: {np.min(df_solps['R'].values)}")
                xpt_spar_slist.append(float(df_solps["Spar"][local_Rmin_solps]))
    
                for case_name, case_obj in cs.items():
                    print(f"Plotting {case_name} with {title}...")
                    # case_name = list(cs.keys())[0] # if input multiple files, select the first one
                    ds = cs[case_name].ds.isel(t=-1)
                    df = get_1d_poloidal_data(ds, params=[p[0] for p in items], region=region_pol, sepadd=ring, target_first=True, interpolate_midplane=False, interpolate_radial=False)
        
                    # Find Hermes-3's x-point
                    local_Rmin_h = np.argmin(df["R"].values) # return location
                    # print(f"Hermes location of Rmin: {local_Rmin_h}, value: {np.min(df['R'].values)}")
                    xpt_spar_hlist.append(float(df["Spar"][local_Rmin_h]))
                    # print(f"The coordinate is: {df["Spar"][local_Rmin_h]}")
            
                    # Hermes-3 
                    axi.plot(np.abs(df["Spar"]), np.abs(df[param]), label=f"Hermes-3 {case_name}")
                    axi.set_title(title)
                    axi.set_ylabel(ylabel)
                    axi.set_xlabel("")
                    axi.grid(True, alpha=0.7)
                    
                    if logy:
                        axi.set_yscale("log")
    
            xpt_h = min(xpt_spar_hlist)
            xpt_s = min(xpt_spar_hlist)
            # print(xpt_h)
           
            # Plot x-point vertical line
            for i, axi in enumerate(ax.flat):
                axi.axvline(xpt_h, color='r', alpha = 0.2)
                axi.axvline(xpt_s, color='r', alpha = 0.2, ls = '--')
            
            # Legend setting
            format_legend_and_axes(fig, ax, nitems, "$S_{\\parallel}$")
            fig.savefig(f"{figures_png_path}/benchmark_rings_{category}.png")



def plot_bm_profiles_radial(cs, bal, region_rad, figures_png_path):
    plots = {
        "state_var": [
            ("Te",  "e temperature", "T [eV]",            False),
            ("Td+", "d+ temperature","T [eV]",            False),
            ("Td",  "d temperature", "T [eV]",            False),
            ("Ne",  "e density",     "density [$m^{-3}$]", True),
            ("Nd+", "d+ density",    "density [$m^{-3}$]", True),
            ("Nd",  "d density",     "density [$m^{-3}$]", True),
            ("Pe",  "e pressure",    "pressure [Pa]",      True),
            ("Pd+", "d+ pressure",   "pressure [Pa]",      True),
            ("Pd",  "d pressure",    "pressure [Pa]",      True),
        ],
        "reactions": [
            ("Sd+_rec", "Recombination", "Rate [$m^{-3} s^{-1}$]", True),
            ("Sd+_iz",  "Ionisation",           "Rate [$m^{-3} s^{-1}$]",          True),
            ("Edd+_cx", "Charge exchange",      "Energy [$W m^{3}$]",   True),
            # ("Fdd+_cx", "Charge exchange",      "Momentum [$kg m^{-2} s^{-2}$]",   True),
            # ("NVd",     "d Parallel momentum",  "Momentum [$kg m^{-2} s^{-2}$]",   True),
            ("NVd+",    "d+ Parallel momentum", "Momentum [$kg m^{-2} s^{-2}$]",   True),
        ]
    }
   
    for category, items in plots.items():
        nitems = len(items)
        print(f"Number of plots = {nitems} in category {category}")
        
        cols = 3
        rows = math.ceil(nitems / cols)
        fig, ax = plt.subplots(rows, cols, figsize=(12, 4 * rows), squeeze=False)

        # Loop through all the profiles
        for idx, (param, title, ylabel, logy) in enumerate(items):

            r, c = divmod(idx, 3)
            axi = ax[r,c]

            # Read SOLPS first (only once if multi Hermes-3 files)
            # SOLPS dist = Srad
            df_solps = bal.get_1d_radial_data(params=[p[0] for p in items], region=region_rad)
            axi.plot(df_solps["dist"], np.abs(df_solps[param]), label=f"SOLPS", ls = "--")

            # Loop through all profiles for each file
            for case_name, case_obj in cs.items():
                print(f"Plotting {case_name} with {title}...")
                ds = cs[case_name].ds.isel(t=-1)
                df = get_1d_radial_data(ds, params=[p[0] for p in items], region=region_rad)

                # Hermes-3 
                axi.plot(df["Srad"], np.abs(df[param]), label=f"Hermes-3 {case_name}")
        
                axi.set_title(title)
                axi.set_ylabel(ylabel)
                axi.set_xlabel("")
                axi.grid(True, alpha=0.7)
            
                if logy:
                    axi.set_yscale("log")

        format_legend_and_axes(fig, ax, nitems, "$X-X_{sep}$ [m]")
        fig.savefig(f"{figures_png_path}/benchmark_radial_{category}.png")

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
        {"vars": ["Ne", "Nd+", "Nd"], "ylabel": "Density (m^-3)", "yscale": "linear", "title_prefix": ""},
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
    figures_png_path = args.output
    if not os.path.exists(f"./{figures_png_path}"):
        os.makedirs(figures_png_path)

    sepadd_array = [1]
    plot_bm_profiles_radial(case, balance, args.region_rad, figures_png_path)
    plot_bm_profiles_fieldline(case, balance, args.region_pol, sepadd_array, figures_png_path)
    # plot_bm_plasma_overlap(case, balance, args.region_pol, args.region_rad, figures_png_path)

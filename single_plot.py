
import matplotlib
import pandas as pd
import numpy as np
import matplotlib

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


#%load_ext autoreload
#%autoreload 2
print("Done")

def read_file():

    # Read file 

    db = CaseDB(
        case_dir = r"/users/jpm590/2dspace/run/",
        grid_dir = r"/users/jpm590/2dspace/hermes-3/build-mc-master"
    )
    
    toload = [
        dict(name="MAST-U", id="251007-2D-MASTU", unnormalise_geom = True, use_xhermes = True, squash = True)
    ]
    cs = {}
    for case in toload:
        cs[case["name"]] = db.load_case_2D(case["id"], use_squash = case["squash"], verbose = True)
        cs[case["name"]].extract_2d_tokamak_geometry()
    m = cs["MAST-U"].ds.metadata
    print(f'Species in model: \n {m["species"]}')
    print(f'\nCharged species: \n {m["charged_species"]}')

    return cs

def make_plot_srad(cs):

    ds = cs["MAST-U"].ds.isel(t=-1)
    fig, ax = plt.subplots(figsize = (4,4))

    df_midplane = get_1d_radial_data(ds, params = ["Te", "Td+", "Td"], region = "omp")
    ax.plot(df_midplane["Srad"], df_midplane["Te"], label = "Te")
    ax.plot(df_midplane["Srad"], df_midplane["Td+"], label = "Td+")
    ax.plot(df_midplane["Srad"], df_midplane["Td"], label = "Td")
    ax.set_xlabel("$X-X_{sep}$ [m]")
    ax.set_ylabel("Temperature [eV]")
    ax.set_title("Midplane temperature")
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/midplane_temperature.png")

    fig, ax = plt.subplots(figsize = (4,4))
    df_midplane = get_1d_radial_data(ds, params = ["Ne", "Nd+", "Nd"], region = "omp")
    ax.plot(df_midplane["Srad"], df_midplane["Ne"], label = "Ne")
    ax.plot(df_midplane["Srad"], df_midplane["Nd+"], label = "Nd+")
    ax.plot(df_midplane["Srad"], df_midplane["Nd"], label = "Nd")
    ax.set_xlabel("$X-X_{sep}$ [m]")
    ax.set_ylabel("density [eV]")
    ax.set_title("Midplane density")
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/midplane_density.png")

def make_plot_spar(cs):

    '''
    
    Plot the whole parallel distance
    function: get_1d_poloidal_data
    Must specify which "ring" of the fieldline by
    1. sepadd: offset from separatrix
    2. sepdist: distance from separartrix 

    '''
        
    ds = cs["MAST-U"].ds.isel(t=-1)
    fig, ax = plt.subplots(figsize = (4,4))


    df_fieldline = get_1d_poloidal_data(ds, params = ["Te", "Td+", "Td"], region = "outer_lower", sepdist =  0.001)
    ax.plot(df_fieldline["Spar"], df_fieldline["Te"], label = "Te")
    ax.plot(df_fieldline["Spar"], df_fieldline["Td+"], label = "Td+")
    ax.plot(df_fieldline["Spar"], df_fieldline["Td"], label = "Td")
    ax.set_xlabel("$S_{\\parallel}$ [m]")
    ax.set_ylabel("Temperature [eV]")
    ax.set_title("Temperatures")
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/fieldline_temperature.png")

    fig, ax = plt.subplots(figsize = (4,4))
    df_fieldline = get_1d_poloidal_data(ds, params = ["Ne", "Nd+", "Nd"], region = "outer_lower", sepdist =  0.001)
    ax.plot(df_fieldline["Spar"], df_fieldline["Ne"], label = "Ne")
    ax.plot(df_fieldline["Spar"], df_fieldline["Nd+"], label = "Nd+")
    ax.plot(df_fieldline["Spar"], df_fieldline["Nd"], label = "Nd")
    ax.set_xlabel("$S_{\\parallel}$ [m]")
    ax.set_ylabel("density [1/m3]")
    ax.set_title("densitys")
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/fieldline_density.png")


def make_plot(cs):
    
    coord_list = ["Spar", "Srad"]
    
    # Group parameters by type
    param_groups = {
        "Temperature": ["Te", "Td", "Td+"],
        "Density": ["Ne", "Nd", "Nd+"],
        "Pressure": ["Pe", "Pd", "Pd+"],
        "Ionisation": ["Sd+_iz"],
        "Recombination": ["Sd+_rec"]
        "Charge_exchange": ["Fdd+_cx"]
        "Impurity": ["Sd_pump"]
    }
    
    ds = cs["MAST-U"].ds.isel(t=-1)
    
    for coord in coord_list:
        for group_name, param_list in param_groups.items():
            print(f"Plotting ... {group_name}")
            fig, ax = plt.subplots(figsize=(4, 4))

            # Get the right dataframe depending on coordinate
            if coord == "Spar":
                df = get_1d_poloidal_data(ds, params=param_list, region="outer_lower", sepdist=0.001)
                x_label = "$S_{\\parallel}$ [m]"
                x_key = "Spar"
            else:
                if group_name == "Ionisation" or group_name == "Charge_exchange" or group_name == "Recombination":
                    continue
                df = get_1d_radial_data(ds, params=param_list, region="omp")
                x_label = "$X-X_{sep}$ [m]"
                x_key = "Srad"
    
            # Plot all parameters in this group
            for param in param_list:
                ax.plot(df[x_key], df[param], label=param)
    
            # Axis and labels
            ax.set_xlabel(x_label)
            if group_name == "Temperature":
                ax.set_ylabel("Temperature [eV]")
            elif group_name == "Density":
                ax.set_ylabel("Density [m$^{-3}$]")
            elif group_name == "Pressure":
                ax.set_ylabel("Pressure [Pa]")
            elif group_name == "Ionisation":
                ax.set_ylabel("Density transfer rate [$m^{-3}s^{-1}$]")
                ax.set_title("Ionisation")
            elif group_name == "Recombination":
                ax.set_ylabel("Density transfer rate [$m^{-3}s^{-1}$]")
                ax.set_title("Recombination")
            elif group_name == "Charge_exchange":
                ax.set_ylabel("Momentum transfer rate [$kg \cdot m^{-2}s^{-2}$]")
                ax.set_title("Charge exchange")
    
            ax.set_title(f"{group_name.capitalize()}s vs {coord}")
            ax.legend()
            fig.tight_layout()
    
            # Save with descriptive name
            fig.savefig(f"figures/{coord}_{group_name}.png")
            plt.close(fig)


def main():
    case = read_file()
    # make_plot_srad(case)
    # make_plot_spar(case)
    make_plot(case)

if __name__ == "__main__":
    
    start_time = time.time()
    main()
    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds")
    

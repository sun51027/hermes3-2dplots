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

figures_pdf_path = args.output + "_figures_pdf"
figures_png_path = args.output + "_figures_png"

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

    m = cs["MAST-U"].ds.metadata
    # print(f'Species in model: \n {m["species"]}')
    # print(f'\nCharged species: \n {m["charged_species"]}')

    ds = cs["MAST-U"].ds.isel(t=-1)
    # print(list(ds.data_vars))
    # print(ds)

    ''' All variable list'''
    # print("All variable list")
    # for name, da in ds.data_vars.items():
    #     print(f"{name}")

    '''Show specific variable data'''
    np.set_printoptions(threshold=np.inf)
    name = "Rc"
    # print(f"Specific variable's data: {name}")
    # print(ds[name])
    # print(ds[name].values)
    # print(ds[name].load())

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
    fig.savefig(f"{figures_pdf_path}/midplane_temperature.pdf")
    fig.savefig(f"{figures_png_path}/midplane_temperature.png")

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
    fig.savefig(f"{figures_pdf_path}/midplane_density.pdf")
    fig.savefig(f"{figures_png_path}/midplane_density.png")

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


    df_fieldline = get_1d_poloidal_data(ds, params = ["Te", "Td+", "Td"], region = "outer_lower", sepadd = 5)
    ax.plot(df_fieldline["Spar"].values[::-1], df_fieldline["Te"], label = "Te")
    ax.plot(df_fieldline["Spar"].values[::-1], df_fieldline["Td+"], label = "Td+")
    ax.plot(df_fieldline["Spar"].values[::-1], df_fieldline["Td"], label = "Td")
    ax.set_xlabel("$S_{\\parallel}$ [m]")
    ax.set_ylabel("Temperature [eV]")
    ax.set_title("Temperatures")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"figures_pdf/fieldline_temperature.pdf")
    # fig.savefig(f"figures_png/fieldline_temperature.png")

    # fig, ax = plt.subplots(figsize = (4,4))
    # df_fieldline = get_1d_poloidal_data(ds, params = ["Ne", "Nd+", "Nd"], region = "outer_lower", sepdist =  0.001)
    # ax.plot(df_fieldline["Spar"], df_fieldline["Ne"], label = "Ne")
    # ax.plot(df_fieldline["Spar"], df_fieldline["Nd+"], label = "Nd+")
    # ax.plot(df_fieldline["Spar"], df_fieldline["Nd"], label = "Nd")
    # ax.set_xlabel("$S_{\\parallel}$ [m]")
    # ax.set_ylabel("density [1/m3]")
    # ax.set_title("densitys")
    # ax.legend()
    # fig.tight_layout()
    # fig.savefig(f"{figures_pdf_path}/fieldline_density.pdf")
    # fig.savefig(f"{figures_png_path}/fieldline_density.png")


def make_plot(cs, region_rad, region_pol, idx_ring):
   
    # The X-axis is from target in LHS to upstream RHS

    coord_list = ["Spar", "Srad"]
    
    # Group parameters by type
    param_groups = {
        "Temperature": ["Te", "Td", "Td+"],
        "Density": ["Ne", "Nd", "Nd+"],
        "Pressure": ["Pe", "Pd", "Pd+"],
        "Ionisation": ["Sd+_iz"],
        "Recombination": ["Sd+_rec"],
        "Charge_exchange": ["Fdd+_cx"],
        # "Impurity": ["Sd_pump"],
        # "Radiation_cooling":["Rc"]
        
    }
    
    ds = cs["MAST-U"].ds.isel(t=-1)
    
    for coord in coord_list:
        for group_name, param_list in param_groups.items():
            print(f"Plotting ... {coord} {group_name}")
            fig, ax = plt.subplots(figsize=(4, 4))

            # Get the right dataframe depending on coordinate
            if coord == "Spar":
                # Fieldline
                df = get_1d_poloidal_data(ds, params=param_list, region=region_pol , sepadd=idx_ring)
                x_label = "$S_{\\parallel}$ [m]"
                x_key = "Spar"
                x_name = "fieldline"
            else:
                # if group_name == "Ionisation" or group_name == "Charge_exchange" or group_name == "Recombination":
                #     continue
                # Radial profile
                df = get_1d_radial_data(ds, params=param_list, region=region_rad)
                # df = get_1d_radial_data(ds, params=param_list, region="omp")
                x_label = "$X-X_{sep}$ [m]"
                x_key = "Srad"
                x_name = "radial"
    
            # Plot all parameters in this group
            for param in param_list:
                if coord == "Spar":
                    # Fieldline 
                    ax.plot(df[x_key].values[::-1], np.abs(df[param]), label=param)
                    if args.scale == "log":
                        ax.set_yscale("log")
                    else:
                        ax.set_yscale("linear")
                else: 
                    # Radial profile
                    ax.plot(df[x_key], df[param], label=param)
                    if args.scale == "log":
                        ax.set_yscale("log")
                    else:
                        ax.set_yscale("linear")
    
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
                ax.set_ylabel("Momentum transfer rate [$kg \\cdot m^{-2}s^{-2}$]")
                ax.set_title("Charge exchange")
    
            ax.legend()
            ax.grid(True, alpha=0.5)
            fig.tight_layout()
    
            # Save with descriptive name
            fig.savefig(f"{figures_png_path}/{x_name}_{group_name}.png")
            fig.savefig(f"{figures_pdf_path}/{x_name}_{group_name}.pdf")
            plt.close(fig)

def make_plot_diff_coeff(cs):

    ds = cs["MAST-U"].ds.isel(t=-1)
    fig, ax = plt.subplots(figsize = (4,4))
    df_fieldline = get_1d_poloidal_data(ds, 
            params = ["anomalous_Chi_d+","anomalous_Chi_e", "anomalous_D_d+", "anomalous_D_e", "anomalous_nu_d+", "anomalous_nu_e"], 
            region = "outer_lower", sepdist =  0.001)

    print(f"Chi_d+ = {df_fieldline["anomalous_Chi_d+"]}")
    print(f"Chi_e  = {df_fieldline["anomalous_Chi_e"]}")
    print(f"D_d+   = {df_fieldline["anomalous_D_d+"]}")
    print(f"D_e    = {df_fieldline["anomalous_D_e"]}")
    print(f"nu_d+  = {df_fieldline["anomalous_nu_d+"]}")
    print(f"nu_e   = {df_fieldline["anomalous_nu_e"]}")
    ax.plot(df_fieldline["Spar"], df_fieldline["anomalous_Chi_d+"], label = "$\\chi_{d+}$")
    ax.plot(df_fieldline["Spar"], df_fieldline["anomalous_Chi_e"], label = "$\\chi_{e}$")
    ax.plot(df_fieldline["Spar"], df_fieldline["anomalous_D_d+"], label = "$D_{d+}$")
    ax.plot(df_fieldline["Spar"], df_fieldline["anomalous_D_e"], label = "$D_{e}$")
    ax.plot(df_fieldline["Spar"], df_fieldline["anomalous_nu_d+"], label = "$\\nu_{d+}$")
    ax.plot(df_fieldline["Spar"], df_fieldline["anomalous_nu_e"], label = "$\\nu_{e}$")
    ax.set_xlabel("$S_{\\parallel}$ [m]")
    ax.set_ylabel("Anomalous coefficient [m/s]")
    ax.set_title("")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{figures_png_path}/fieldline_diff_coeff.png")
    fig.savefig(f"{figures_pdf_path}/fieldline_diff_coeff.pdf")


    df_midplane = get_1d_radial_data(ds, 
            params = ["anomalous_Chi_d+","anomalous_Chi_e", "anomalous_D_d+", "anomalous_D_e", "anomalous_nu_d+", "anomalous_nu_e"], 
            region = "omp")

    print(f"Chi_d+ = {df_midplane["anomalous_Chi_d+"]}")
    print(f"Chi_e  = {df_midplane["anomalous_Chi_e"]}")
    print(f"D_d+   = {df_midplane["anomalous_D_d+"]}")
    print(f"D_e    = {df_midplane["anomalous_D_e"]}")
    print(f"nu_d+  = {df_midplane["anomalous_nu_d+"]}")
    print(f"nu_e   = {df_midplane["anomalous_nu_e"]}")
    ax.plot(df_midplane["Spar"], df_midplane["anomalous_Chi_d+"], label = "$\\chi_{d+}$")
    ax.plot(df_midplane["Spar"], df_midplane["anomalous_Chi_e"], label = "$\\chi_{e}$")
    ax.plot(df_midplane["Spar"], df_midplane["anomalous_D_d+"], label = "$D_{d+}$")
    ax.plot(df_midplane["Spar"], df_midplane["anomalous_D_e"], label = "$D_{e}$")
    ax.plot(df_midplane["Spar"], df_midplane["anomalous_nu_d+"], label = "$\\nu_{d+}$")
    ax.plot(df_midplane["Spar"], df_midplane["anomalous_nu_e"], label = "$\\nu_{e}$")
    ax.set_xlabel("$S_{\\parallel}$ [m]")
    ax.set_ylabel("Anomalous coefficient [m/s]")
    ax.set_title("")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{figures_png_path}/midplane_diff_coeff.png")
    fig.savefig(f"{figures_pdf_path}/midplane_diff_coeff.pdf")

def plot_Lz_function(cs):
    ds = cs["MAST-U"].ds.isel(t=-1)
    fig, ax = plt.subplots(1,2, figsize=(10,4))
    Lz = abs(ds["Rc"].hermesm.clean_guards())/(ds["Ne"].hermesm.clean_guards()*ds["Ne"].hermesm.clean_guards()*0.02)
    
    ax[0].scatter(ds["Te"], Lz, s=5)
    ax[0].set_xlim(0,40)
    # ax[0].set_ylim(0,9e-32)
    ax[0].set_title("The whole domain")
    ax[0].grid(True, alpha =0.5)
    ax[0].set_yscale("log")
    ax[0].set_ylabel("L_z function [W/m^3]")
    ax[0].set_xlabel("$T_e$ [eV]")

    df_fieldline = get_1d_poloidal_data(ds, params = ["Rc", "Ne", "Te","Spar"], region = "outer_lower", sepadd = 2)

    Lz = abs(df_fieldline["Rc"])/(df_fieldline["Ne"]*df_fieldline["Ne"]*0.02)
    ax[1].set_title("2nd SOL ring at outer lower region")
    ax[1].scatter(df_fieldline["Te"], Lz, s=5)
    ax[1].set_xlim(0,40)
    ax[1].set_yscale("log")

    ax[1].grid(True, alpha =0.5)
    ax[1].set_xlabel("$T_e$ [eV]")

    fig.savefig(f"{figures_pdf_path}/Lz_function.pdf")

# def plot_polygon(cs):


# ------------------ main ---------------------

def main():

    print(args.scale)
    case = read_file(args.input)

    '''For test'''
    # make_plot_spar(case)

    '''Plot all function'''
    # region_rad = "omp"
    # region_pol = "outer_lower"
    # idx_ring = 15 # deps on nx, e.g. 0 - 19 
    # idx_ring = args.sepadd
    make_plot(case, args.region_rad, args.region_pol, args.sepadd)

    # make_plot_diff_coeff(case)
    # plot_Lz_function(case)

if __name__ == "__main__":
    
    start_time = time.time()
    if not os.path.exists(f"./{figures_pdf_path}"):
        os.makedirs(figures_pdf_path)
    if not os.path.exists(f"./{figures_png_path}"):
        os.makedirs(figures_png_path)
    main()
    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds")
    

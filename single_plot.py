
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
        dict(name="250929-shorter-test", id="250929-shorter-test", unnormalise_geom = True, use_xhermes = True, squash = True)
    ]
    cs = {}
    for case in toload:
        cs[case["name"]] = db.load_case_2D(case["id"], use_squash = case["squash"], verbose = True)
        cs[case["name"]].extract_2d_tokamak_geometry()
    m = cs["250929-shorter-test"].ds.metadata
    print(f'Species in model: \n {m["species"]}')
    print(f'\nCharged species: \n {m["charged_species"]}')

    return cs

def make_plot_srad(cs):

    ds = cs["250929-shorter-test"].ds.isel(t=-1)
    fig, ax = plt.subplots(figsize = (4,4))
    # plot_selection(ds, ds.hermesm.select_region("outer_upper_target"))
    # plt.savefig('test.png')
    # plot_selection(ds, ds.hermesm.select_region("outer_lower"))
    # plt.savefig('core.png')
    df_midplane = get_1d_radial_data(ds, params = ["Te", "Td+", "Td"], region = "omp")
    ax.plot(df_midplane["Srad"], df_midplane["Te"], label = "Te")
    ax.plot(df_midplane["Srad"], df_midplane["Td+"], label = "Td+")
    ax.plot(df_midplane["Srad"], df_midplane["Td"], label = "Td")
    ax.set_xlabel("$X-X_{sep}$ [m]")
    ax.set_ylabel("Temperature [eV]")
    ax.set_title("Midplane temperatures")
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/midplane_temperature.png")


def make_plot_spar(cs):

    '''
    
    Plot the whole parallel distance
    function: get_1d_poloidal_data
    Must specify which "ring" of the fieldline by
    1. sepadd: offset from separatrix
    2. sepdist: distance from separartrix 

    '''

    coord = ["Spar", "Srad"]
    temp = ["Te", "Td", "Td+"]
    
    for i in coord:
        
    ds = cs["250929-shorter-test"].ds.isel(t=-1)
    fig, ax = plt.subplots(figsize = (4,4))
        if coord == "Spar":
            df_fieldline = get_1d_poloidal_data(ds, params = temp, region = "outer_lower", sepdist =  0.001)
        else:
            df_midplane = get_1d_radial_data(ds, params = temp, , region = "omp")
    ax.plot(df_fieldline["Spar"], df_fieldline["Te"], label = "Te")
    ax.plot(df_fieldline["Spar"], df_fieldline["Td+"], label = "Td+")
    ax.plot(df_fieldline["Spar"], df_fieldline["Td"], label = "Td")
    ax.set_xlabel("$S_{\parallel}$ [m]")
    ax.set_ylabel("Temperature [eV]")
    ax.set_title("Temperatures")
    ax.legend()
    fig.tight_layout()
    fig.savefig("figures/fieldline_temperature.png")

def main():
    case = read_file()
    make_plot_srad(case)
    make_plot_spar(case)

if __name__ == "__main__":
    main()


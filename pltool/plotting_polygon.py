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

def plot_temperature(cs, figures_png_path):

    for case_name, case_obj in cs.items():
        print(f"Plotting T polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)
        fig, ax = plt.subplots(1,3 , figsize=(10,6), constrained_layout=True)

        ds["Te"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True,) #vmin=1, vmax=100)
        ds["Td"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True,) #vmin=1, vmax=100)
        ds["Td+"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True,)# vmin=1, vmax=100)

        fig.suptitle(f"{case_name}")
        fig.savefig(f"{figures_png_path}/{case_name}_T_polygon.png")

def plot_density(cs, figures_png_path):

    for case_name, case_obj in cs.items():
        print(f"Plotting N polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)
        fig, ax = plt.subplots(1,3 , figsize=(10,6), constrained_layout=True)

        ds["Ne"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
        ds["Nd"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
        ds["Nd+"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True,)
        fig.suptitle(f"{case_name}")
        fig.savefig(f"{figures_png_path}/{case_name}_N_polygon.png")

def plot_pressure(cs, figures_png_path):

    for case_name, case_obj in cs.items():
        print(f"Plotting P polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)
        fig, ax = plt.subplots(1,3 , figsize=(10,6), constrained_layout=True)

        ds["Pe"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
        ds["Pd"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
        ds["Pd+"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True,)
        fig.suptitle(f"{case_name}")
        fig.savefig(f"{figures_png_path}/{case_name}_P_polygon.png")

def plot_momentum(cs, figures_png_path):

    for case_name, case_obj in cs.items():
        print(f"Plotting P polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)
        fig, ax = plt.subplots(1,3 , figsize=(10,6), constrained_layout=True)

        ds["NVe"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
        ds["NVd"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
        ds["NVd+"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True,)
        fig.suptitle(f"{case_name}")
        fig.savefig(f"{figures_png_path}/{case_name}_NV_polygon.png")
        
def plot_momentum(cs, figures_png_path):

    for case_name, case_obj in cs.items():
        print(f"Plotting P polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)
        fig, ax = plt.subplots(1,3 , figsize=(10,6), constrained_layout=True)

        ds["NVe"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
        ds["NVd"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
        ds["NVd+"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True,)
        fig.suptitle(f"{case_name}")
        fig.savefig(f"{figures_png_path}/{case_name}_NV_polygon.png")
def run_polygon():

    setup_matplotlib()
    parser = build_base_parser()
    args = parser.parse_args()
    case= read_files(args.input)

    ## create output directory
    figures_png_path = args.output + "_figures_png"
    if not os.path.exists(f"./{figures_png_path}"):
        os.makedirs(figures_png_path)


    ## Run
    plot_density(case, figures_png_path)
    plot_temperature(case, figures_png_path)
    plot_pressure(case, figures_png_path)

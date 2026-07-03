from .common import build_base_parser, read_files, setup_matplotlib
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

def _plot_solps_polygon(bal, solps_param, ax, logscale=False):
    """Plot a single SOLPS field on the given axis, matching the Hermes-3 style.

    solps_param is None when the field has no SOLPS counterpart; in that case
    the axis is left blank so the row still lines up with the Hermes-3 column.
    """
    if solps_param is None or bal is None:
        ax.set_axis_off()
        return
    bal.plot_2d(solps_param, ax=ax, cmap="Spectral_r", antialias=True,
                logscale=logscale, separatrix=True)


def plot_temperature(cs, figures_png_path, solps=False, bal=None):

    # (hermes_param, solps_param); solps_param=None -> no SOLPS counterpart
    fields = [
        ("Te",  "Te"),
        ("Td",  "Td"),
        ("Td+", "Td+"),
    ]

    for case_name, case_obj in cs.items():
        print(f"Plotting T polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        if solps:
            fig, ax = plt.subplots(2, len(fields), figsize=(4 * len(fields), 8),
                                   constrained_layout=True, squeeze=False)
            for col, (h_param, s_param) in enumerate(fields):
                ds[h_param].bout.polygon(ax=ax[0, col], cmap="Spectral_r", antialias=True)
                _plot_solps_polygon(bal, s_param, ax[1, col])
            ax[0, 0].set_ylabel("Hermes-3")
            ax[1, 0].set_ylabel("SOLPS")
        else:
            fig, ax = plt.subplots(1, 3, figsize=(10, 6), constrained_layout=True)
            ds["Te"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True,) #vmin=1, vmax=100)
            ds["Td"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True,) #vmin=1, vmax=100)
            ds["Td+"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True,)# vmin=1, vmax=100)

        fig.suptitle(f"{case_name}")
        fig.savefig(f"{figures_png_path}/{case_name}_T_polygon.png")

def plot_density(cs, figures_png_path, solps=False, bal=None):

    # SOLPS: Nd -> Na (already aliased as "Nd"); Nd+ -> Ne (quasi-neutrality)
    fields = [
        ("Ne",  "Ne"),
        ("Nd",  "Nd"),
        ("Nd+", "Ne"),
    ]

    for case_name, case_obj in cs.items():
        print(f"Plotting N polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        if solps:
            fig, ax = plt.subplots(2, len(fields), figsize=(4 * len(fields), 8),
                                   constrained_layout=True, squeeze=False)
            for col, (h_param, s_param) in enumerate(fields):
                ds[h_param].bout.polygon(ax=ax[0, col], cmap="Spectral_r", antialias=True, logscale=True)
                _plot_solps_polygon(bal, s_param, ax[1, col], logscale=True)
            ax[0, 0].set_ylabel("Hermes-3")
            ax[1, 0].set_ylabel("SOLPS")
        else:
            fig, ax = plt.subplots(1, 3, figsize=(10, 6), constrained_layout=True)
            ds["Ne"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
            ds["Nd"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
            ds["Nd+"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True,)

        fig.suptitle(f"{case_name}")
        fig.savefig(f"{figures_png_path}/{case_name}_N_polygon.png")

def plot_pressure(cs, figures_png_path, solps=False, bal=None):

    # SOLPS: Pd -> Pa (atom pressure)
    fields = [
        ("Pe",  "Pe"),
        ("Pd",  "Pa"),
        ("Pd+", "Pd+"),
    ]

    for case_name, case_obj in cs.items():
        print(f"Plotting P polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        if solps:
            fig, ax = plt.subplots(2, len(fields), figsize=(4 * len(fields), 8),
                                   constrained_layout=True, squeeze=False)
            for col, (h_param, s_param) in enumerate(fields):
                ds[h_param].bout.polygon(ax=ax[0, col], cmap="Spectral_r", antialias=True, logscale=True)
                _plot_solps_polygon(bal, s_param, ax[1, col], logscale=True)
            ax[0, 0].set_ylabel("Hermes-3")
            ax[1, 0].set_ylabel("SOLPS")
        else:
            fig, ax = plt.subplots(1, 3, figsize=(10, 6), constrained_layout=True)
            ds["Pe"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
            ds["Pd"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
            ds["Pd+"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True,)

        fig.suptitle(f"{case_name}")
        fig.savefig(f"{figures_png_path}/{case_name}_P_polygon.png")

def plot_momentum(cs, figures_png_path, solps=False, bal=None):

    # SOLPS only provides d+ parallel momentum; NVe/NVd have no counterpart
    fields = [
        ("NVe",  None),
        ("NVd",  None),
        ("NVd+", "NVd+"),
    ]

    for case_name, case_obj in cs.items():
        print(f"Plotting NV polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        if solps:
            fig, ax = plt.subplots(2, len(fields), figsize=(4 * len(fields), 8),
                                   constrained_layout=True, squeeze=False)
            for col, (h_param, s_param) in enumerate(fields):
                ds[h_param].bout.polygon(ax=ax[0, col], cmap="Spectral_r", antialias=True, logscale=True, vmax=1e-3, vmin=-1e-3)
                _plot_solps_polygon(bal, s_param, ax[1, col], logscale=True)
            ax[0, 0].set_ylabel("Hermes-3")
            ax[1, 0].set_ylabel("SOLPS")
        else: 
            fig, ax = plt.subplots(1, 3, figsize=(10, 6), constrained_layout=True)
            ds["NVe"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
            ds["NVd"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
            ds["NVd+"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True,)

        fig.suptitle(f"{case_name}") 
        fig.savefig(f"{figures_png_path}/{case_name}_NV_polygon.png")
# def run_polygon():

#     setup_matplotlib()
#     parser = build_base_parser()
#     args = parser.parse_args()
#     case= read_files(args.input)

#     ## create output directory
#     figures_png_path = args.output + "_figures_png"
#     if not os.path.exists(f"./{figures_png_path}"):
#         os.makedirs(figures_png_path)


#     ## Run
#     plot_density(case, figures_png_path)
#     plot_temperature(case, figures_png_path)
#     plot_pressure(case, figures_png_path)

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

def _as_solps_dict(bal):
    """Normalise the SOLPS argument to {name: SOLPScase}.

    read_solps_balance() returns a dict (one entry per case, like read_files),
    but a bare SOLPScase or None is also accepted for convenience.
    """
    if bal is None:
        return {}
    if isinstance(bal, dict):
        return bal
    return {"SOLPS": bal}


def _solps_grid(nfields, bals):
    """Figure with one Hermes-3 row on top and one row per SOLPS case below."""
    nrows = 1 + len(bals)
    fig, ax = plt.subplots(nrows, nfields, figsize=(4 * nfields, 4 * nrows),
                           constrained_layout=True, squeeze=False)
    return fig, ax


def _label_rows(ax, bals):
    ax[0, 0].set_ylabel("Hermes-3")
    for row, sname in enumerate(bals, start=1):
        ax[row, 0].set_ylabel(sname)


def _plot_solps_polygon(bal, solps_param, ax, logscale=False, cmap="Spectral_r",
                        vmin=None, vmax=None):
    """Plot a single SOLPS field on the given axis, matching the Hermes-3 style.

    solps_param is None when the field has no SOLPS counterpart; in that case
    the axis is left blank so the row still lines up with the Hermes-3 column.
    """
    if solps_param is None or bal is None:
        ax.set_axis_off()
        return
    bal.plot_2d(solps_param, ax=ax, cmap=cmap, antialias=True,
                logscale=logscale, vmin=vmin, vmax=vmax, separatrix=True)


def _log_vlims(values, vmin, vmax):
    """Default (vmin, vmax) for a log plot of a (possibly signed) field.

    Uses the data's natural range, but keeps the bounds off exactly zero so
    xbout's norm is valid regardless of sign:
      * non-negative data -> both bounds positive  (LogNorm)
      * mixed-sign data    -> vmin < 0 < vmax       (SymLogNorm, negatives shown)
      * non-positive data  -> vmin < 0, tiny +vmax  (SymLogNorm straddling zero)
    A zero bound would make SymLogNorm's linthresh zero ("linthresh must be
    positive"); a signed bound would make LogNorm invalid. Explicit vmin/vmax
    pass through unchanged so they can be adjusted freely.
    """
    vals = np.asarray(values)
    finite = vals[np.isfinite(vals)]
    lo = float(finite.min()) if finite.size else 0.0
    hi = float(finite.max()) if finite.size else 0.0
    nonzero = np.abs(finite[finite != 0])
    floor = float(nonzero.min()) if nonzero.size else 1.0  # smallest magnitude present
    if vmin is None:
        vmin = floor if lo >= 0 else lo
    if vmax is None:
        vmax = hi if hi > 0 else floor
    return vmin, vmax


def _plot_hermes_reaction(da, ax, vmin=None, vmax=None):
    """Log-plot a Hermes-3 reaction field over its natural range (vmin/vmax optional)."""
    da = da.hermes.clear_guards()
    vmin, vmax = _log_vlims(da.values, vmin, vmax)
    da.bout.polygon(ax=ax, cmap="Spectral_r", antialias=True,
                    logscale=True, vmin=vmin, vmax=vmax)


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
            bals = _as_solps_dict(bal)
            fig, ax = _solps_grid(len(fields), bals)
            for col, (h_param, s_param) in enumerate(fields):
                ds[h_param].hermes.clear_guards().bout.polygon(ax=ax[0, col], cmap="Spectral_r", antialias=True)
                for row, b in enumerate(bals.values(), start=1):
                    _plot_solps_polygon(b, s_param, ax[row, col], logscale=False)
            _label_rows(ax, bals)
        else:
            fig, ax = plt.subplots(1, 3, figsize=(10, 6), constrained_layout=True)
            ds["Te"].hermes.clear_guards().bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True,) #vmin=1, vmax=100)
            ds["Td"].hermes.clear_guards().bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True,) #vmin=1, vmax=100)
            ds["Td+"].hermes.clear_guards().bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True,)# vmin=1, vmax=100)

        fig.suptitle(f"{case_name}")
        fig.savefig(f"{figures_png_path}/{case_name}_T_polygon.png")

def plot_density(cs, figures_png_path, solps=False, bal=None):

    # SOLPS: Nd -> Na (already aliased as "Nd"); Nd+ -> Ne (quasi-neutrality)
    fields = [
        ("Ne",  "Ne"),
        ("Nd",  "Nd"),
        ("Nd+", "Nd+"),
    ]

    for case_name, case_obj in cs.items():
        print(f"Plotting N polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        if solps:
            bals = _as_solps_dict(bal)
            fig, ax = _solps_grid(len(fields), bals)
            for col, (h_param, s_param) in enumerate(fields):
                ds[h_param].hermes.clear_guards().bout.polygon(ax=ax[0, col], cmap="Spectral_r", antialias=True, logscale=False)
                for row, b in enumerate(bals.values(), start=1):
                    _plot_solps_polygon(b, s_param, ax[row, col], logscale=False )
            _label_rows(ax, bals)
        else:
            fig, ax = plt.subplots(1, 3, figsize=(10, 6), constrained_layout=True)
            ds["Ne"].hermes.clear_guards().bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = False, vmin=1e18 , vmax=3e19)
            ds["Nd"].hermes.clear_guards().bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = False, vmin=1e18 , vmax=3e19)
            ds["Nd+"].hermes.clear_guards().bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = False,vmin=1e18 , vmax=3e19)

        fig.suptitle(f"{case_name}")
        fig.savefig(f"{figures_png_path}/{case_name}_N_polygon.png")

def plot_pressure(cs, figures_png_path, solps=False, bal=None):

    # SOLPS: Pd -> Pa (atom pressure)
    fields = [
        ("Pe",  "Pe"),
        ("Pd",  "Pd"),
        ("Pd+", "Pd+"),
    ]

    for case_name, case_obj in cs.items():
        print(f"Plotting P polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        if solps:
            bals = _as_solps_dict(bal)
            fig, ax = _solps_grid(len(fields), bals)
            for col, (h_param, s_param) in enumerate(fields):
                ds[h_param].hermes.clear_guards().bout.polygon(ax=ax[0, col], cmap="Spectral_r", antialias=True, logscale=True)
                for row, b in enumerate(bals.values(), start=1):
                    _plot_solps_polygon(b, s_param, ax[row, col], logscale=True)
            _label_rows(ax, bals)
        else:
            fig, ax = plt.subplots(1, 3, figsize=(10, 6), constrained_layout=True)
            ds["Pe"].hermes.clear_guards().bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
            ds["Pd"].hermes.clear_guards().bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
            ds["Pd+"].hermes.clear_guards().bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True,)

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
            bals = _as_solps_dict(bal)
            fig, ax = _solps_grid(len(fields), bals)
            for col, (h_param, s_param) in enumerate(fields):
                ds[h_param].hermes.clear_guards().bout.polygon(ax=ax[0, col], cmap="Spectral_r", antialias=True, logscale=True, vmax=1e-3, vmin=-1e-3)
                for row, b in enumerate(bals.values(), start=1):
                    _plot_solps_polygon(b, s_param, ax[row, col], logscale=True)
            _label_rows(ax, bals)
        else: 
            fig, ax = plt.subplots(1, 3, figsize=(10, 6), constrained_layout=True)
            ds["NVe"].hermes.clear_guards().bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
            ds["NVd"].hermes.clear_guards().bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
            ds["NVd+"].hermes.clear_guards().bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True,)

        fig.suptitle(f"{case_name}") 
        fig.savefig(f"{figures_png_path}/{case_name}_NV_polygon.png")

def plot_reaction(cs, figures_png_path, solps=False, bal=None, vmin=None, vmax=None):

    # SOLPS only provides d+ parallel momentum; NVe/NVd have no counterpart
    fields = [
        ("Sd+_rec", "Sd+_rec"),
        ("Sd+_iz",  "Sd+_iz"),
        ("Edd+_cx", "Edd+_cx"),
    ]

    for case_name, case_obj in cs.items():
        print(f"Plotting reaction polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        if solps:
            bals = _as_solps_dict(bal)
            fig, ax = _solps_grid(len(fields), bals)
            for col, (h_param, s_param) in enumerate(fields):
                _plot_hermes_reaction(ds[h_param], ax[0, col], vmin=vmin, vmax=vmax)
                for row, b in enumerate(bals.values(), start=1):
                    if s_param is not None:
                        s_vmin, s_vmax = _log_vlims(b.bal[s_param], vmin, vmax)
                    else:
                        s_vmin, s_vmax = vmin, vmax
                    _plot_solps_polygon(b, s_param, ax[row, col], logscale=True,
                                        vmin=s_vmin, vmax=s_vmax)
            _label_rows(ax, bals)
        else:
            fig, ax = plt.subplots(1, 3, figsize=(10, 6), constrained_layout=True)
            for col, (h_param, _) in enumerate(fields):
                _plot_hermes_reaction(ds[h_param], ax[col], vmin=vmin, vmax=vmax)

        fig.suptitle(f"{case_name}")
        fig.savefig(f"{figures_png_path}/{case_name}_reaction_polygon.png")
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

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


def _limits(vmin, vmax, same_limit):
    """Colour limits for one field.


    same_limit=True  -> the (vmin, vmax) written in the fields table, so the
                        Hermes-3 row and every SOLPS row share one colour scale
                        and can actually be compared by eye.
    same_limit=False -> (None, None), i.e. every axis autoscales as before.
    """
    return (vmin, vmax) if same_limit else (None, None)


def _solps_grid(nfields, bals):
    """Figure with one Hermes-3 row on top and one row per SOLPS case below."""
    nrows = 1 + len(bals)
    fig, ax = plt.subplots(nrows, nfields, figsize=(2.5 * nfields, 4 * nrows),
                           constrained_layout=True, squeeze=False)
    return fig, ax


def _label_rows(ax, bals):
    ax[0, 0].set_ylabel("Hermes-3")
    for row, sname in enumerate(bals, start=1):
        ax[row, 0].set_ylabel(sname)


# Font sizes for the polygon figures only; the 1D profile plots keep the
# global rcParams set in setup_matplotlib().
POLYGON_FONT = dict(title=12, label=11, tick=9, suptitle=15)


def _apply_fonts(fig, font=None):
    """Resize every text element of a polygon figure.

    Neither bout.polygon() nor SOLPScase.plot_2d() takes a fontsize argument,
    so the sizes are applied afterwards, once the axes (including the colorbar
    axes plot_2d appends) already exist. Pass e.g. font=dict(tick=7) to
    override individual entries of POLYGON_FONT.
    """
    f = dict(POLYGON_FONT, **(font or {}))
    for axi in fig.axes:
        axi.title.set_fontsize(f["title"])
        axi.xaxis.label.set_fontsize(f["label"])
        axi.yaxis.label.set_fontsize(f["label"])
        axi.tick_params(axis="both", which="both", labelsize=f["tick"])
    if fig._suptitle is not None:
        fig._suptitle.set_fontsize(f["suptitle"])


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
                logscale=logscale, vmin=vmin, vmax=vmax, separatrix=True, )


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


def plot_temperature(cs, figures_png_path, solps=False, bal=None, same_limit=False, font=None):

    # (hermes_param, solps_param, vmin, vmax); solps_param=None -> no SOLPS counterpart
    # vmin/vmax are only used when same_limit=True; tune them per machine/case.
    fields = [
        ("Td",  "Td",  0, 500),
        ("Te",  "Te",  0, 500),
        ("Td+", "Td+", 0, 500),
    ]

    bals = _as_solps_dict(bal) if solps else {}

    for case_name, case_obj in cs.items():
        print(f"Plotting T polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        fig, ax = _solps_grid(len(fields), bals)
        for col, (h_param, s_param, fvmin, fvmax) in enumerate(fields):
            vmin, vmax = _limits(fvmin, fvmax, same_limit)
            ds[h_param].hermes.clear_guards().bout.polygon(
                ax=ax[0, col], cmap="Spectral_r", antialias=True,
                logscale=False, vmin=vmin, vmax=vmax)
            for row, b in enumerate(bals.values(), start=1):
                _plot_solps_polygon(b, s_param, ax[row, col], logscale=False,
                                    vmin=vmin, vmax=vmax)
        _label_rows(ax, bals)

        fig.suptitle(f"{case_name}")
        _apply_fonts(fig, font)
        fig.savefig(f"{figures_png_path}/{case_name}_T_polygon.png")

def plot_density(cs, figures_png_path, solps=False, bal=None, same_limit=False, font=None):

    # SOLPS: Nd -> Na (already aliased as "Nd"); Nd+ -> Ne (quasi-neutrality)
    fields = [
        ("Nd",  "Nd",  1e11, 5e20),
        ("Ne",  "Ne",  1e17, 1e19),
        ("Nd+", "Nd+", 1e17, 1e19),
    ]

    bals = _as_solps_dict(bal) if solps else {}

    for case_name, case_obj in cs.items():
        print(f"Plotting N polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        fig, ax = _solps_grid(len(fields), bals)
        for col, (h_param, s_param, fvmin, fvmax) in enumerate(fields):
            vmin, vmax = _limits(fvmin, fvmax, same_limit)
            ds[h_param].hermes.clear_guards().bout.polygon(
                ax=ax[0, col], cmap="Spectral_r", antialias=True,
                logscale=True, vmin=vmin, vmax=vmax)
            for row, b in enumerate(bals.values(), start=1):
                _plot_solps_polygon(b, s_param, ax[row, col], logscale=True,
                                    vmin=vmin, vmax=vmax)
        _label_rows(ax, bals)

        fig.suptitle(f"{case_name}")
        _apply_fonts(fig, font)
        fig.savefig(f"{figures_png_path}/{case_name}_N_polygon.png")

def plot_pressure(cs, figures_png_path, solps=False, bal=None, same_limit=False, font=None):

    # SOLPS: Pd -> Pa (atom pressure)
    fields = [
        ("Pd",  "Pd",  1e-6, 1e1),
        ("Pe",  "Pe",  1e-2, 1e2),
        ("Pd+", "Pd+", 1e-2, 1e2),
    ]

    bals = _as_solps_dict(bal) if solps else {}

    for case_name, case_obj in cs.items():
        print(f"Plotting P polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        fig, ax = _solps_grid(len(fields), bals)
        for col, (h_param, s_param, fvmin, fvmax) in enumerate(fields):
            vmin, vmax = _limits(fvmin, fvmax, same_limit)
            ds[h_param].hermes.clear_guards().bout.polygon(
                ax=ax[0, col], cmap="Spectral_r", antialias=True,
                logscale=True, vmin=vmin, vmax=vmax)
            for row, b in enumerate(bals.values(), start=1):
                _plot_solps_polygon(b, s_param, ax[row, col], logscale=True,
                                    vmin=vmin, vmax=vmax)
        _label_rows(ax, bals)

        fig.suptitle(f"{case_name}")
        _apply_fonts(fig, font)
        fig.savefig(f"{figures_png_path}/{case_name}_P_polygon.png")

def plot_momentum(cs, figures_png_path, solps=False, bal=None, same_limit=False, font=None):

    # SOLPS only provides d+ parallel momentum; NVe/NVd have no counterpart.
    # Momentum changes sign, so the limits straddle zero (SymLogNorm).
    fields = [
        ("NVe",  None,   -1e-3, 1e-3),
        # ("NVd+", "NVd+", -1e-3, 1e-3),
        ("NVd",  None,   -1e-3, 1e-3),
    ]

    bals = _as_solps_dict(bal) if solps else {}

    for case_name, case_obj in cs.items():
        print(f"Plotting NV polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        fig, ax = _solps_grid(len(fields), bals)
        for col, (h_param, s_param, fvmin, fvmax) in enumerate(fields):
            vmin, vmax = _limits(fvmin, fvmax, same_limit)
            ds[h_param].hermes.clear_guards().bout.polygon(
                ax=ax[0, col], cmap="Spectral_r", antialias=True,
                logscale=True, vmin=vmin, vmax=vmax)
            for row, b in enumerate(bals.values(), start=1):
                _plot_solps_polygon(b, s_param, ax[row, col], logscale=True,
                                    vmin=vmin, vmax=vmax)
        _label_rows(ax, bals)

        fig.suptitle(f"{case_name}")
        _apply_fonts(fig, font)
        fig.savefig(f"{figures_png_path}/{case_name}_NV_polygon.png")

def plot_reaction(cs, figures_png_path, solps=False, bal=None, vmin=None, vmax=None,
                  same_limit=False, font=None):

    # Reaction rates span many decades, hence the log limits below
    fields = [
        ("Sd+_rec", "Sd+_rec", -1e18, 1e18),
        ("Sd+_iz",  "Sd+_iz",  1e16, 1e24),
        # ("Edd+_cx", "Edd+_cx", 1e2,  1e6),
    ]

    bals = _as_solps_dict(bal) if solps else {}

    for case_name, case_obj in cs.items():
        print(f"Plotting reaction polygon for {case_name}")

        ds = case_obj.ds.isel(t=-1)

        fig, ax = _solps_grid(len(fields), bals)
        for col, (h_param, s_param, fvmin, fvmax) in enumerate(fields):
            cvmin, cvmax = _limits(fvmin, fvmax, same_limit)
            # explicit vmin/vmax arguments still override the table
            cvmin = vmin if vmin is not None else cvmin
            cvmax = vmax if vmax is not None else cvmax

            _plot_hermes_reaction(ds[h_param], ax[0, col], vmin=cvmin, vmax=cvmax)
            for row, b in enumerate(bals.values(), start=1):
                if s_param is not None and not same_limit:
                    # autoscale each SOLPS panel over its own range
                    s_vmin, s_vmax = _log_vlims(b.bal[s_param], cvmin, cvmax)
                else:
                    s_vmin, s_vmax = cvmin, cvmax
                _plot_solps_polygon(b, s_param, ax[row, col], logscale=True,
                                    vmin=s_vmin, vmax=s_vmax)
        _label_rows(ax, bals)

        fig.suptitle(f"{case_name}")
        _apply_fonts(fig, font)
        fig.savefig(f"{figures_png_path}/{case_name}_reaction_polygon.png")

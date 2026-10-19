from .common import (
    build_base_parser, read_files, setup_matplotlib, read_solps_balance,
    format_legend_and_axes, H3_COLORS, H3_STYLES, SOLPS_COLORS, SOLPS_STYLES,
)
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



def core_avg_hermes(ds, param, region="core_noguards"):
    """Volume-weighted mean of `param` over closed flux surfaces (inside separatrix)."""
    sub = ds.hermesm.select_region(region)   # radial 0..ixseps1, core poloidal range
    q, dv = sub[param], sub["dv"]
    return float((q * dv).sum() / dv.sum())


def core_avg_solps(bal, param):
    """Volume-weighted mean of `param` inside the separatrix (SOLPS balance)."""
    data = np.asarray(bal.bal[param])        # (poloidal, radial)
    vol = np.asarray(bal.bal["vol"])
    pol = bal.psel["core"]                    # core poloidal indices (inner + outer)
    rad = slice(1, bal.g["sep"])              # inside separatrix, drop inner guard ring
    q = data[pol][:, rad]
    w = vol[pol][:, rad]
    return float(np.nansum(q * w) / np.nansum(w))

# def target_flux_hermes(cs, param, region="outer_lower_target_guards"):
#     ds = cs[case_name].ds.isel(t=-1)
#     sub = ds.hermesm.select_region(region)   # radial 0..ixseps1, core poloidal range
#     dv = sub["dv"]
#     print(f"dv shape = {dv.shape}")
    # df = get_1d_radial_data(ds, params=[p[0] for p in items], region=region_rad)


    

def print_core_averages(cs, bals, params=("Te", "Td+", "Td","Ne", "Nd+","Nd"), region="core_noguards"):
    """Print a volume-weighted core (inside-separatrix) average table: SOLPS vs Hermes-3.

    `bals` is the {name: SOLPScase} dict from read_solps_balance(); a bare
    SOLPScase is still accepted. Percentages are relative to the FIRST SOLPS
    case, which acts as the reference.
    """
    # Accept both a dict of cases and a single SOLPScase
    if not isinstance(bals, dict):
        bals = {"SOLPS": bals}

    ref_name = next(iter(bals)) if bals else None
    print("\n=== Volume-weighted core averages (inside separatrix) ===")
    if ref_name is not None:
        print(f"(percentages are relative to SOLPS {ref_name})")

    header = f"{'param':8s}"
    for sname in bals:
        header += f" | {'SOLPS ' + sname:>20s}"
    for case_name in cs:
        header += f" | {'H3 ' + case_name:>20s}"
    print(header)

    for param in params:
        s_vals = {}
        for sname, bal in bals.items():
            try:
                s_vals[sname] = core_avg_solps(bal, param)
            except Exception as e:
                s_vals[sname] = float("nan")
                print(f"  [warn] SOLPS {sname} {param}: {e}")

        ref = s_vals.get(ref_name, float("nan"))

        def cell(val):
            rel = (val - ref) / ref * 100.0 if ref else float("nan")
            return f" | {val:12.4g} ({rel:+.1f}%)"

        row = f"{param:8s}"
        for sname in bals:
            row += cell(s_vals[sname])

        for case_name in cs:
            ds = cs[case_name].ds.isel(t=-1)
            try:
                h_val = core_avg_hermes(ds, param, region=region)
            except Exception as e:
                h_val = float("nan")
                print(f"  [warn] Hermes-3 {case_name} {param}: {e}")
            row += cell(h_val)
        print(row)
    print("=" * 60)



def xpt_spar(df, absolute=False):
    """Parallel coordinate of the x-point, i.e. Spar where R along the ring is smallest."""
    imin = int(np.argmin(df["R"].values))
    spar = df["Spar"].values[imin]
    return float(abs(spar) if absolute else spar)


def ring_kwargs(ds, region_pol, sepadd=None, sepdist=None,
                target_first=True, interpolate_midplane=False):
    """Put both codes on the SAME flux surface.

    Hermes-3 uses its own ring index; SOLPS is radially interpolated to that
    ring's distance from the separatrix, which is how interpolateSOLPStoHermes
    matched them in the first place. Give either sepadd or sepdist.
    """
    radial_region = "omp" if "outer" in region_pol else "imp"
    rad = get_1d_radial_data(ds, params=["R"], region=radial_region, guards=False)
    sol = rad[rad["region"] == "sol"]

    if (sepadd is None) == (sepdist is None):
        raise ValueError("Give either sepadd or sepdist, not both")

    if sepdist is not None:
        idx = (sol["Srad"] - sepdist).abs().idxmin()
        sepadd = sol.index.get_loc(idx)

    sepdist_used = sol.iloc[sepadd]["Srad"]
    print(f"[ring] sepadd = {sepadd}, sepdist = {sepdist_used * 1e3:.2f} mm")

    hermes_kw = dict(region=region_pol, sepadd=sepadd, target_first=target_first,
                     interpolate_midplane=interpolate_midplane,
                     interpolate_radial=False)
    solps_kw = dict(region=f"{region_pol}_extra" if interpolate_midplane else region_pol,
                    sepdist=sepdist_used, target_first=target_first,
                    interpolate_midplane=interpolate_midplane,
                    interpolate_radial=False, extrapolate_radial=False)
    return hermes_kw, solps_kw, sepadd, sepdist_used


def plot_bm_profiles_fieldline(cs, bals, region_pol, idx_ring_array, figures_png_path,
                               ring_mode="sepadd"):

    # param, title, ylabel, logy, symlog
    plots = {
        "state_var": [
            ("Te",  "e temperature", "T [eV]",            True,False),
            ("Ne",  "e density",     "density [$m^{-3}$]", True,False),
            ("Pe",  "e pressure",    "pressure [Pa]",      True,False),
            ("Td+", "d+ temperature","T [eV]",            False,False),
            ("Nd+", "d+ density",    "density [$m^{-3}$]", True,False),
            ("Pd+", "d+ pressure",   "pressure [Pa]",      True,False),
            ("Td",  "d temperature", "T [eV]",            False,False),
            ("Nd",  "d density",     "density [$m^{-3}$]", True,False),
            ("Pd",  "d pressure",    "pressure [Pa]",      True,False),
        ],
        "momentum":[
            ("Vd",      "d Parallel velocity",    "Velocity [$m s^{-1}$]",       False  , True),
            ("NVd",     "d Parallel momentum",  "Momentum [$kg m^{-2} s^{-1}$]", False  , True),
            # ("NVd+",    "d+ Parallel momentum", "Momentum [$kg m^{-2} s^{-1}$]", False  , True),
        ],
        "reactions": [
            ("Sd+_rec", "Recombination", "Source [$m^{-3} s^{-1}$]",        True, True),
            ("Sd+_iz",  "Ionisation",           "Source [$m^{-3} s^{-1}$]", True, False),
            # ("Edd+_cx", "Ion energy source due to charge exchange",  "Energy [$W m^{-3}$]",    True, False),
            # ("Ed+_iz",  "Ion energy source due to ionisation ",      "Energy [$W m^{-3}$]",    True, False),
        ]
    }

    for ring in idx_ring_array:

        for category, items in plots.items():
            params = [p[0] for p in items]
            nitems = len(items)
            print(f"Number of plots = {nitems} in category {category}")
            cols = 3
            rows = math.ceil(nitems / cols)
            fig, ax = plt.subplots(rows, cols, figsize=(10, 3.5 * rows), squeeze=False)

            # Read each case once per ring/category, not once per parameter.
            # Both codes are matched onto the same flux surface, see ring_kwargs.
            ref_ds = cs[next(iter(cs))].ds.isel(t=-1)
            h_kw, s_kw, sepadd, sepdist_used = ring_kwargs(ref_ds, region_pol,
                                                           **{ring_mode: ring})

            solps_dfs = {n: b.get_1d_poloidal_data(params=params, **s_kw) for n, b in bals.items()}

            h3_dfs = {}
            for case_name, case_obj in cs.items():
                print(f"Reading {case_name} ring {sepadd}...")
                h3_dfs[case_name] = get_1d_poloidal_data(case_obj.ds.isel(t=-1), params=params, **h_kw)

            # X-point = Spar where R is smallest, on each code's own x-axis
            # xpt_s = xpt_spar(df_solps)
            xpt_s = min(xpt_spar(df) for df in solps_dfs.values())
            xpt_h = min(xpt_spar(df, absolute=True) for df in h3_dfs.values())

            for idx, (param, title, ylabel, logy, symlog) in enumerate(items):
                r, c = divmod(idx, cols)
                axi = ax[r,c]

                # SOLPS
                for j, (sname, df) in enumerate(solps_dfs.items()):
                    # axi.plot(df["Spol"], (df[param]), label=f"SOLPS {sname}",
                    axi.plot(df["Spar"], (df[param]), label=f"SOLPS {sname}",
                             color=SOLPS_COLORS[j % len(SOLPS_COLORS)],
                             **SOLPS_STYLES[j % len(SOLPS_STYLES)])

                # Hermes-3
                for i, (case_name, df) in enumerate(h3_dfs.items()):
                    # axi.plot((df["Spol"]), (df[param]), label=f"Hermes-3 {case_name}", **H3_STYLES[i % len(H3_STYLES)])
                    if params=="Edd+_cx":
                        axi.plot((df["Spar"]), np.abs(df[param]), label=f"Hermes-3 {case_name}", **H3_STYLES[i % len(H3_STYLES)])
                    else:
                        axi.plot((df["Spar"]), (df[param]), label=f"Hermes-3 {case_name}", **H3_STYLES[i % len(H3_STYLES)])

                # X-point vertical lines
                # axi.axvline(xpt_h, color='r', alpha = 0.2)
                # axi.axvline(xpt_s, color='r', alpha = 0.2, ls = '--')

                axi.set_title(title)
                axi.set_ylabel(ylabel)
                axi.set_xlabel("$S_{\\parallel}$")
                axi.grid(True, alpha=0.7)

                if logy:
                    axi.set_yscale("log")
                    # axi.set_xscale("log")
                if symlog:
                    axi.set_yscale("symlog")


            # Legend setting
            # format_legend_and_axes(fig, ax, nitems, "$S_{pol}$ [m]")
            format_legend_and_axes(fig, ax, nitems, "$S_{\\parallel}$ [m]")
            # fig.suptitle(f"ring {sepadd}  (sepdist = {sepdist_used * 1e3:.2f} mm)")
            fig.savefig(f"{figures_png_path}/benchmark_ring{sepadd}_{category}.png")



def plot_bm_profiles_radial(cs, bals, region_rad, figures_png_path):
    plots = {
        "state_var": [
            ("Te",  "e temperature", "T [eV]",            False),
            ("Ne",  "e density",     "density [$m^{-3}$]", True),
            ("Pe",  "e pressure",    "pressure [Pa]",      True),
            ("Td+", "d+ temperature","T [eV]",            False),
            ("Nd+", "d+ density",    "density [$m^{-3}$]", True),
            ("Pd+", "d+ pressure",   "pressure [Pa]",      True),
            ("Td",  "d temperature", "T [eV]",            False),
            ("Nd",  "d density",     "density [$m^{-3}$]", True),
            ("Pd",  "d pressure",    "pressure [Pa]",      True),
        ],
        "momentum":[
            ("Vd",      "d Parallel velocity",    "Velocity [$m s^{-1}$]",       True),
            ("NVd",     "d Parallel momentum",  "Momentum [$kg m^{-2} s^{-1}$]", True  ),
            # ("NVd+",    "d+ Parallel momentum", "Momentum [$kg m^{-2} s^{-1}$]", False),
        ],
        "reactions": [
            ("Sd+_rec", "Ion density source due to recombination",        "$m^{-3} s^{-1}$", True),
            ("Sd+_iz",  "Ion desnity source due to ionisation",           "$m^{-3} s^{-1}$",          True),
            ("Edd+_cx", "Ion energy source due to charge exchange \n (-Edd+_cx)",  "Energy source [$W m^{-3}$]",   True),
            ("Ed+_iz",  "Ion energy source due to ionisation \n (Ed+_iz) ",      "Energy source [$W m^{-3}$]",    True),
            ("Edd+_cxiz",  "Ion energy source due to iz & cx \n ",      "$W m^{-3}$",    True),
            # ("Edd+_cxiz",  "Ion energy source due to iz & cx \n (-Edd+_cx + Ed+_iz)",      "$W m^{-3}$",    True),
            # ("Edd+_cx", "Neutral energy source due to charge exchange \n (Edd+_cx)",  "Energy source [$W m^{-3}$]",   True),
            # ("Edd+_cxiz",  "Neutral energy source due to iz & cx \n ",      "Energy source[$W m^{-3}$]",    True),
            # ("Edd+_cxiz",  "Neutral energy source due to iz & cx \n (Edd+_cx - Ed+_iz)",      "Energy source[$W m^{-3}$]",    True),
            # ("Fdd+_cx", "Charge exchange",      "Momentum [$kg m^{-2} s^{-2}$]",   True),
        ]
    }
   
    for category, items in plots.items():
        nitems = len(items)
        print(f"Number of plots = {nitems} in category {category}")

        # Variables to load: plotted ones + dependencies of derived quantities
        params = [p[0] for p in items]
        if "Edd+_cxiz" in params:
            params += [p for p in ("Edd+_cx", "Ed+_iz") if p not in params]

        cols = 3
        rows = math.ceil(nitems / cols)
        fig, ax = plt.subplots(rows, cols, figsize=(10, 3.5 * rows), squeeze=False)

        # Loop through all the profiles
        for idx, (param, title, ylabel, logy) in enumerate(items):

            r, c = divmod(idx, 3)
            axi = ax[r,c]

            # SOLPS loop
            for j, (sname, bal) in enumerate(bals.items()):

                df_solps = bal.get_1d_radial_data(params=[p[0] for p in items], region=region_rad)
                
                if param=="Edd+_cx" or param=="Ed+_iz":
                    continue
                elif param=="Edd+_cxiz":
                    # To show "neutral eneryg source",
                    # need to flip the sign of eirene_mc_eapl_shi_bal (ion energy srouce)
                    y_s = -df_solps[param]
                else:
                    y_s = df_solps[param]

                axi.plot(df_solps["dist"], y_s,
                         label=f"SOLPS {sname}",
                         color=SOLPS_COLORS[j % len(SOLPS_COLORS)],
                         **SOLPS_STYLES[j % len(SOLPS_STYLES)])
            # df_solps = bal.get_1d_radial_data(params=[p[0] for p in items], region=region_rad)
            # axi.plot(df_solps["dist"], np.abs(df_solps[param]), label="SOLPS", color=SOLPS_COLOR, **SOLPS_STYLE)

            # Hermes-3 loop
            for i, (case_name, case_obj) in enumerate(cs.items()):
                print(f"Plotting {case_name} with {title}...")
                ds = cs[case_name].ds.isel(t=-1)
                df = get_1d_radial_data(ds, params=params, region=region_rad)

                # Hermes-3
                if param=="Edd+_cx":
                    # Show in "neutral energy source"
                    y = df[param]

                    # If you want to show "ion source"
                    # y = -df[param]
                elif param=="Edd+_cxiz":
                    # Show in "neutral energy source"
                    y = df["Edd+_cx"] - df["Ed+_iz"]

                    # If you want to show "ion source"
                    # y = -df["Edd+_cx"] + df["Ed+_iz"]
                else:
                    y = df[param]

                axi.plot(df["Srad"], y, label=f"Hermes-3 {case_name}", **H3_STYLES[i % len(H3_STYLES)])
        
                axi.set_title(title)
                axi.set_ylabel(ylabel)
                axi.set_xlabel("$X-X_{sep}$ [m]")
                axi.grid(True, alpha=0.7)
            
                if logy:
                    axi.set_yscale("symlog")
                    # axi.set_yscale("log")

        format_legend_and_axes(fig, ax, nitems, "$X-X_{sep}$ [m]")
        fig.savefig(f"{figures_png_path}/benchmark_radial_{category}.png")

def plot_bm_profiles_fieldline_multiring(cs, bals, region_pol, idx_ring_array,
                                         figures_png_path,
                                         alpha_step=0.2, alpha_min=0.2,
                                         normalise_spar=False):
    """Overlay several SOL rings on one axis: colour = code/case, alpha = ring.

    idx_ring_array : list of sepadd values, e.g. [1, 2, 3]
                     first ring gets alpha=1.0, then 0.8, 0.6 ...
    """

    plots = {
        "state_var": [
            ("Te",  "e temperature", "T [eV]",             False),
            ("Td+", "d+ temperature","T [eV]",             False),
            ("Td",  "d temperature", "T [eV]",             False),
            ("Ne",  "e density",     "density [$m^{-3}$]",  True),
            ("Nd+", "d+ density",    "density [$m^{-3}$]",  True),
            ("Nd",  "d density",     "density [$m^{-3}$]",  True),
            ("Pe",  "e pressure",    "pressure [Pa]",       True),
            ("Pd+", "d+ pressure",   "pressure [Pa]",       True),
            ("Pd",  "d pressure",    "pressure [Pa]",       True),
        ],
        "reactions": [
            ("Sd+_rec", "Recombination",        "Rate [$m^{-3} s^{-1}$]",        True),
            ("Sd+_iz",  "Ionisation",           "Rate [$m^{-3} s^{-1}$]",        True),
            ("Edd+_cx", "Charge exchange",      "Energy [$W m^{3}$]",            True),
            ("NVd+",    "d+ Parallel momentum", "Momentum [$kg m^{-2} s^{-1}$]", True),
        ],
    }

    rings = list(idx_ring_array)
    # 1.0, 0.8, 0.6 ... never fading below alpha_min
    ring_alpha = {r: max(1.0 - alpha_step * i, alpha_min) for i, r in enumerate(rings)}

    for category, items in plots.items():
        params = [p[0] for p in items]
        nitems = len(items)
        cols = 3
        rows = math.ceil(nitems / cols)
        fig, ax = plt.subplots(rows, cols, figsize=(12, 3.5 * rows), squeeze=False)

        # ---- read every (case, ring) once ------------------------------------
        solps_dfs, h3_dfs = {}, {}
        for ring in rings:
            kw = dict(region=region_pol, sepadd=ring, target_first=True,
                      interpolate_midplane=False, interpolate_radial=False)

            for sname, b in bals.items():
                try:
                    solps_dfs[(sname, ring)] = b.get_1d_poloidal_data(params=params, **kw)
                except Exception as e:
                    print(f"  [warn] SOLPS {sname} ring {ring}: {e}")

            for case_name, case_obj in cs.items():
                print(f"Reading {case_name} ring {ring}...")
                try:
                    h3_dfs[(case_name, ring)] = get_1d_poloidal_data(
                        case_obj.ds.isel(t=-1), params=params, **kw)
                except Exception as e:
                    print(f"  [warn] Hermes-3 {case_name} ring {ring}: {e}")

        def xaxis(df, absolute=False):
            s = np.abs(df["Spar"]) if absolute else df["Spar"]
            return s / np.nanmax(np.abs(s)) if normalise_spar else s

        # ---- plot -------------------------------------------------------------
        for idx, (param, title, ylabel, logy) in enumerate(items):
            r, c = divmod(idx, cols)
            axi = ax[r, c]

            # SOLPS
            for (sname, ring), df in solps_dfs.items():
                j = list(bals).index(sname)
                style = dict(SOLPS_STYLES[j % len(SOLPS_STYLES)])   # copy: it already holds alpha
                style["alpha"] = ring_alpha[ring]
                axi.plot(xaxis(df), (df[param]),
                         color=SOLPS_COLORS[j % len(SOLPS_COLORS)],
                         label="_nolegend_", **style)

            # Hermes-3
            for (case_name, ring), df in h3_dfs.items():
                i = list(cs).index(case_name)
                style = dict(H3_STYLES[i % len(H3_STYLES)])
                style["alpha"] = ring_alpha[ring]
                axi.plot(xaxis(df, absolute=False), df[param],
                         color=H3_COLORS[i % len(H3_COLORS)],
                         label="_nolegend_", **style)

            # X-point of each ring, faded the same way
            for ring in rings:
                hs = [xpt_spar(df, absolute=True)
                      for (cn, rg), df in h3_dfs.items() if rg == ring]
                ss = [xpt_spar(df) for (sn, rg), df in solps_dfs.items() if rg == ring]
                if hs:
                    axi.axvline(min(hs), color="r", alpha=0.2 * ring_alpha[ring])
                if ss:
                    axi.axvline(min(ss), color="r", ls="--", alpha=0.2 * ring_alpha[ring])

            axi.set_title(title)
            axi.set_ylabel(ylabel)
            axi.set_xlabel("")
            axi.grid(True, alpha=0.7)

            if logy:
                axi.set_yscale("log")
                axi.set_xscale("log")

        # ---- two-block legend: one entry per code/case, one per ring ----------
        leg_ax = ax[0, 0]
        for i, case_name in enumerate(cs):
            leg_ax.plot([], [], color=H3_COLORS[i % len(H3_COLORS)],
                        label=f"Hermes-3 {case_name}",
                        **H3_STYLES[i % len(H3_STYLES)])
        for j, sname in enumerate(bals):
            style = dict(SOLPS_STYLES[j % len(SOLPS_STYLES)])
            style["alpha"] = 1.0
            leg_ax.plot([], [], color=SOLPS_COLORS[j % len(SOLPS_COLORS)],
                        label=f"SOLPS {sname}", **style)
        for ring in rings:
            leg_ax.plot([], [], color="k", lw=3, alpha=ring_alpha[ring],
                        label=f"ring {ring} ($\\alpha$={ring_alpha[ring]:.1f})")

        format_legend_and_axes(fig, ax, nitems, "$S_{\\parallel}$")
        tag = "-".join(str(r) for r in rings)
        fig.savefig(f"{figures_png_path}/benchmark_multiring{tag}_{category}.png")


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



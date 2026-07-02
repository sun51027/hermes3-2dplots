from .common import build_base_parser, read_files, setup_matplotlib, read_solps_balance
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

def plot_single_profiles(cs, region_rad, region_pol, idx_ring, args, figures_png_path):
  

    coord_list = ["Spar", "Srad"]
    
    # Group parameters by type
    param_groups = {
        "Temperature": ["Te", "Td", "Td+"],
        "Density": ["Ne", "Nd", "Nd+"],
        "Pressure": ["Pe", "Pd", "Pd+"],
        "Ionisation": ["Sd+_iz"],
        "Recombination": ["Sd+_rec"],
        "Charge_exchange_mom": ["Fdd+_cx"],
        "Charge_exchange_energy": ["Edd+_cx"],
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
            elif group_name == "Charge_exchange_mom":
                ax.set_ylabel("Momentum transfer rate [$kg \\cdot m^{-2}s^{-2}$]")
                ax.set_title("Charge exchange momentum")
            elif group_name == "Charge_exchange_energy":
                ax.set_ylabel("Energy transfer due to CX [$ W/m^{3}$]")
                ax.set_title("Charge exchange energy")
    
            ax.legend()
            ax.grid(True, alpha=0.5)
            fig.tight_layout()
    
            # Save with descriptive name
            fig.savefig(f"{figures_png_path}/{x_name}_{group_name}.png")
#            fig.savefig(f"{figures_pdf_path}/{x_name}_{group_name}.pdf")
            plt.close(fig)

def make_plot_diff_coeff(cs):

    ds = cs["MAST-U"].ds.isel(t=-1)
    fig, ax = plt.subplots(figsize = (4,4))
    df_fieldline = get_1d_poloidal_data(ds, 
            params = ["anomalous_Chi_d+","anomalous_Chi_e", "anomalous_D_d+", "anomalous_D_e", "anomalous_nu_d+", "anomalous_nu_e"], 
            region = "outer_lower", sepdist =  0.001)

    print("\n" + "="*60)
    print(" ANOMALOUS DIFFUSION COEFFICIENTS (Outer Lower Fieldline)")
    print("="*60)
    cols = ["Spar", "anomalous_Chi_d+", "anomalous_Chi_e", "anomalous_D_d+", "anomalous_D_e", "anomalous_nu_d+", "anomalous_nu_e"]
    if all(c in df_fieldline.columns for c in cols):
        print(df_fieldline[cols].to_string(index=False))
    else:
        print(df_fieldline.to_string(index=False))
    print("="*60 + "\n")
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
#    fig.savefig(f"{figures_pdf_path}/fieldline_diff_coeff.pdf")


    df_midplane = get_1d_radial_data(ds, 
            params = ["anomalous_Chi_d+","anomalous_Chi_e", "anomalous_D_d+", "anomalous_D_e", "anomalous_nu_d+", "anomalous_nu_e"], 
            region = "omp")

    print("\n" + "="*60)
    print(" ANOMALOUS DIFFUSION COEFFICIENTS (Outer Midplane)")
    print("="*60)
    cols = ["Spar", "anomalous_Chi_d+", "anomalous_Chi_e", "anomalous_D_d+", "anomalous_D_e", "anomalous_nu_d+", "anomalous_nu_e"]
    if all(c in df_midplane.columns for c in cols):
        print(df_midplane[cols].to_string(index=False))
    else:
        print(df_midplane.to_string(index=False))
    print("="*60 + "\n")
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
#    fig.savefig(f"{figures_pdf_path}/midplane_diff_coeff.pdf")

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

#    fig.savefig(f"{figures_pdf_path}/Lz_function.pdf")

def compare_last_timestep(cs, bal, figures_png_path):
    
    # Dictionary of variables to compare
    # Format: {"name": {"param": "...", "region": "...", "operation": "...", "ylabel": "...", "title": "..."}}
    variables_to_compare = {
        "Te_target_max": {
            "param": "Te",
            "region": "outer_lower_target",
            "operation": "max",
            "ylabel": "$T_{e}^{target,max}$ [eV]",
            "title": "Trend of $T_{e}^{target,max}$ vs puff rate",
            "grad_unit": "eV/s",
            "label_name": "$T_{e}^{target,max}$"
        },
        "Ne_omp_sep": {
            "param": "Ne",
            "region": "omp",
            "operation": "sep",
            "ylabel": "$N_{e}^{omp,sep}$ [$m^{-3}$]",
            "title": "Trend of $N_{e}^{omp,sep}$ vs puff rate",
            "grad_unit": "$m^{-3}s^{-1}$",
            "label_name": "$N_{e}^{omp,sep}$"
        }
    }
    
    keys = list(cs.keys())
    
    # Check if keys can be cast to float 
    # (Change initial value to False if you want to force equally-spaced categorical X-axis)
    x_numeric = False
    for k in keys:
        try:
            float(k)
        except ValueError:
            x_numeric = False
            break
            
    for var_key, var_info in variables_to_compare.items():
        param = var_info["param"]
        region = var_info["region"]
        operation = var_info["operation"]
        ylabel = var_info["ylabel"]
        title = var_info["title"]
        label_name = var_info["label_name"]
        grad_unit = var_info.get("grad_unit", "")
        
        # Get SOLPS value
        df_solps = bal.get_1d_radial_data(params=[param], region=region)
        if operation == "max":
            solps_val = np.abs(df_solps[param]).max()
        elif operation == "sep":
            if "dist" in df_solps.columns:
                idx_sep = np.argmin(np.abs(df_solps["dist"]))
            else:
                idx_sep = np.argmin(np.abs(df_solps["Srad"]))
            solps_val = np.abs(df_solps[param]).iloc[idx_sep]
        else:
            solps_val = 0.0
            
        x_vals = []
        y_vals = []
        gradients = []
        te_firsts = []
        
        fig, ax = plt.subplots(figsize=(6, 6))
        
        for i, (case_name, case_obj) in enumerate(cs.items()):
            ds = case_obj.ds
            
            def get_hermes_val(ds_time):
                df = get_1d_radial_data(ds_time, params=[param], region=region)
                if operation == "max":
                    return np.abs(df[param]).max()
                elif operation == "sep":
                    idx_sep = np.argmin(np.abs(df["Srad"]))
                    return np.abs(df[param]).iloc[idx_sep]
                return 0.0
            
            # Last time step
            ds_last = ds.isel(t=-1)
            val_last = get_hermes_val(ds_last)
            
            # First time step
            ds_first = ds.isel(t=0)
            val_first = get_hermes_val(ds_first)
            
            t_last = float(ds_last["t"].values)
            t_first = float(ds_first["t"].values)
            dt = t_last - t_first
            
            if dt > 0:
                gradient = (val_last - val_first) / dt
            else:
                gradient = 0.0
                
            x_val = float(case_name) if x_numeric else i
            x_vals.append(x_val)
            y_vals.append(val_last)
            gradients.append(gradient)
            te_firsts.append(val_first)
            
        ax.scatter(x_vals, y_vals, marker='o', label=f"Hermes-3 {label_name} (last step)", color='blue')
        
        # Plot horizontal line for SOLPS
        ax.axhline(solps_val, color='r', linestyle='--', label=f"SOLPS {label_name}")
        
        # Add arrows
        for x, y_last, y_first, grad in zip(x_vals, y_vals, te_firsts, gradients):
            if y_last != y_first:
                # Draw arrow from first timestep to last timestep to represent the change
                ax.annotate("", xy=(x, y_last), xytext=(x, y_first),
                            arrowprops=dict(arrowstyle="->", color="black", lw=1.5))
                # Annotate the rate of change (gradient) next to the final point
                # ax.text(x, y_last, f" {grad:.2e} {grad_unit}", ha='left', va='center', fontsize=9)
                
        if not x_numeric:
            ax.set_xticks(range(len(keys)))
            ax.set_xticklabels(keys)
            ax.set_xlabel("Cases")
        else:
            ax.set_xlabel("Puff rate")
            if min(x_vals) > 0 and max(x_vals) / min(x_vals) >= 10:
                ax.set_xscale("log")
                
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.5)
        
        fig.tight_layout()
        fig.savefig(f"{figures_png_path}/compare_last_timestep_{var_key}.png")
        plt.close(fig)

def plot_particle_balance_time(cs, figures_png_path):
    for case_name, case_obj in cs.items():
        ds = case_obj.ds
        
        fig, ax = plt.subplots(figsize=(8, 6))
        # Plasma conservation
        # ddt n_p + n_n = - particleflow_p - particleflow_n  + S_puff - S_pump + S_recycle  
        RHS_sum = 0
        
        print(f"\n{'='*60}")
        print(f" PARTICLE BALANCE REPORT: {case_name}")
        print(f"{'='*60}")
        
        # 1. ddt(total density integral of n+plasma)
        if "Ne" in ds and "Nd" in ds and "dv" in ds:
            # Keep as xarray DataArrays so .differentiate("t") works!
            # Note: We use Ne (electrons) as a proxy for ions (Nd+) due to quasi-neutrality.
            # Adding Ne + Nd tracks the total number of ATOMS.
            total_plasma = (ds["Ne"].hermesm.clean_guards() * ds["dv"]).sum(["x", "theta"])
            total_neutral = (ds["Nd"].hermesm.clean_guards() * ds["dv"]).sum(["x", "theta"])
            total_density_integral = total_plasma + total_neutral
            
            ddt_plasma = total_plasma.differentiate("t")
            ddt_neutral = total_neutral.differentiate("t")
            ddt_total = total_density_integral.differentiate("t")
            
            print("\n[1] Particle Inventory & Time Derivative")
            print(f"    Total plasma (Np)            : {total_plasma.values[-1]:11.4e}")
            print(f"    Total neutral (Nn)           : {total_neutral.values[-1]:11.4e}")
            print(f"    d/dt Np                      : {ddt_plasma.values[-1]:11.4e}")
            print(f"    d/dt Nn                      : {ddt_neutral.values[-1]:11.4e}")
            print(f"    d/dt (Np + Nn)               : {ddt_total.values[-1]:11.4e}")
            
            ax.plot(ds["t"], ddt_total, label="d/dt (Total Particles)", lw=2, color="black", marker="o")
        else:
            print(f"\n[{case_name}] Missing variables for total particle calculation")
            ddt_plasma = ds["t"] * 0
            ddt_neutral = ds["t"] * 0
            ddt_total = ds["t"] * 0
            
        # 2. Puff rate vs time
        if "Sd_src" in ds and "dv" in ds:
            puff_rate = (ds["Sd_src"].hermesm.clean_guards() * ds["dv"]).sum(["x", "theta"])
            ax.plot(ds["t"], puff_rate, label="Puff Rate (Sd_src)", linestyle="--", lw=2)
            
            print("\n[2] Sources & Sinks (1/s)")
            print(f"    Puff rate (Sd_src)           : {puff_rate.values[-1]:11.4e}")
            RHS_sum = RHS_sum + puff_rate
        else:
            puff_rate = ds["t"] * 0
            print("\n[2] Sources & Sinks (1/s)")
            print("    Puff rate (Sd_src)           : Missing")
        
        # 3. Pump rate vs time
        if "Sd_pump" in ds and "dv" in ds:
            pump_rate = np.abs((ds["Sd_pump"].hermesm.clean_guards() * ds["dv"]).sum(["x", "theta"]))
            # Pump is a sink (negative), so we can plot its true value or abs
            ax.plot(ds["t"], np.abs(pump_rate), label="Pump Rate (Sd_pump)", linestyle="--", lw=2)
            print(f"    Pump rate (Sd_pump)          : {pump_rate.values[-1]:11.4e}")
            RHS_sum -= pump_rate
        else:
            pump_rate = ds["t"] * 0
            print("    Pump rate (Sd_pump)          : Missing")
            
        # 4. Target Recycle vs time
        if "Sd_target_recycle" in ds and "dv" in ds:
            target_recycle = (ds["Sd_target_recycle"].hermesm.clean_guards() * ds["dv"]).sum(["x", "theta"])
            # ax.plot(ds["t"], target_recycle, label="Target Recycle", linestyle="--", lw=2)
            
            # If Recycling R = 0.99, then the recycled neutrals = 0.99 * Ion_Outflux
            # The particles permanently LOST to the wall = 0.01 * Ion_Outflux
            # Therefore, Loss = (0.01 / 0.99) * target_recycle
            R = 0.99
            target_loss = ((1 - R) / R) * target_recycle
            ax.plot(ds["t"], target_loss, label="Target Loss = (1-R)/R * Sd_target_recyle, \n R = recycling factor", linestyle="--", lw=2)
            
            print("\n[3] Target Recycling (1/s)")
            print(f"    Target recycle (99% back)    : {target_recycle.values[-1]:11.4e}")
            print(f"    Target loss (1% ion loss)    : {target_loss.values[-1]:11.4e}")
            RHS_sum += target_recycle
        else:
            target_recycle = ds["t"] * 0
            target_loss = ds["t"] * 0
            print("\n[3] Target Recycling (1/s)")
            print("    Target recycle / loss        : Missing")
        
        # 5. Divergence of particle flow
        # particle flow (1/s)= integral volume of div dot (particle flux), so don't need np.gradient 
        # 1/s = m^3 * (1/m) * 1/(m^2 * s), so don't need to integrate dv
        # this means that the flow term is already integrated over the volume
        if "pfd+_tot_xlow" in ds and "pfd+_tot_ylow" in ds and "pfd_adv_par_ylow" in ds and "pfd_adv_perp_xlow" in ds and "pfd_adv_perp_ylow" in ds:
            plasmaflow_x        = (ds["pfd+_tot_xlow"].shift(x=-1) - ds["pfd+_tot_xlow"]).hermesm.clean_guards().sum(["x", "theta"])
            plasmaflow_y        = (ds["pfd+_tot_ylow"].shift(theta=-1) - ds["pfd+_tot_ylow"]).hermesm.clean_guards().sum(["x", "theta"])
            neutralflow_par_y   = (ds["pfd_adv_par_ylow"].shift(theta=-1) - ds["pfd_adv_par_ylow"]).hermesm.clean_guards().sum(["x", "theta"])
            neutralflow_perp_x  = (ds["pfd_adv_perp_xlow"].shift(x=-1) - ds["pfd_adv_perp_xlow"]).hermesm.clean_guards().sum(["x", "theta"])
            neutralflow_perp_y  = (ds["pfd_adv_perp_ylow"].shift(theta=-1) - ds["pfd_adv_perp_ylow"]).hermesm.clean_guards().sum(["x", "theta"])

            total_plasma_flow = plasmaflow_x + plasmaflow_y
            total_neutral_flow = neutralflow_par_y + neutralflow_perp_x + neutralflow_perp_y
            total_flow = total_plasma_flow + total_neutral_flow
            
            print("\n[4] Particle Flows (1/s)")
            print(f"    Total plasma flow            : {total_plasma_flow.values[-1]:11.4e}")
            print(f"    Total neutral flow           : {total_neutral_flow.values[-1]:11.4e}")
            print(f"    Total particle flow          : {total_flow.values[-1]:11.4e}")

            # ax.plot(ds["t"], np.abs(total_flow), label="Total particle flow", linestyle="--", lw=2)
            RHS_sum -= total_flow
        else:
            print("\n[4] Particle Flows (1/s)")
            print("    Particle flow variables      : Missing")
            total_plasma_flow = ds["t"] * 0
            total_neutral_flow = ds["t"] * 0
            total_flow = ds["t"] * 0

        # 6. Recombination rate vs time
        if "Sd+_rec" in ds and "dv" in ds:
            recomb_rate = np.abs((ds["Sd+_rec"].hermesm.clean_guards() * ds["dv"]).sum(["x", "theta"]))
            print("\n[5] Atomic Processes (1/s)")
            print(f"    Recomb rate (Sd+_rec)        : {recomb_rate.values[-1]:11.4e}")
        else:
            recomb_rate = ds["t"] * 0
            print("\n[5] Atomic Processes (1/s)")
            print("    Recomb rate (Sd+_rec)        : Missing")

        # 7. Ionisation rate vs time
        if "Sd+_iz" in ds and "dv" in ds:
            iz_rate = np.abs((ds["Sd+_iz"].hermesm.clean_guards() * ds["dv"]).sum(["x", "theta"]))
            print(f"    Ionisation rate (Sd+_iz)     : {iz_rate.values[-1]:11.4e}")
        else:
            iz_rate = ds["t"] * 0
            print("    Ionisation rate (Sd+_iz)     : Missing")

        # The true particle balance of the entire simulatio        
        # ddt n_p = - particleflow_p + S_iz - S_rec n
        plasma_RHS = -total_plasma_flow + iz_rate - recomb_rate
        plasma_diff = ddt_plasma - plasma_RHS

        # ddt n_n = - particleflow_n - S_iz + S_rec + S_puff - S_pump + S_recycle
        neutral_RHS = -total_neutral_flow - iz_rate + recomb_rate + puff_rate - pump_rate + target_recycle
        neutral_diff = ddt_neutral - neutral_RHS
        
        # ddt n_p + n_n = - particleflow_p - particleflow_n  + S_puff - S_pump + S_recycle  
        # total_RHS = RHS_sum
        total_RHS = plasma_RHS + neutral_RHS
        difference = ddt_total - total_RHS
        
        true_net_balance = puff_rate - pump_rate - target_loss
        
        print("\n[6] Overall Balance (1/s)")
        print("--- Plasma Conservation ---")
        print("    Eq: ddt(N_p) = - Div(Flow_p) + S_iz - S_rec")
        print(f"    - Div(Flow_p)                : {-total_plasma_flow.values[-1]:11.4e}")
        print(f"    + S_iz                       : {iz_rate.values[-1]:11.4e}")
        print(f"    - S_rec                      : {-recomb_rate.values[-1]:11.4e}")
        print(f"    LHS [d/dt(N_p)]              : {ddt_plasma.values[-1]:11.4e}")
        print(f"    RHS [Flows & Sources]        : {plasma_RHS.values[-1]:11.4e}")
        print(f"    Difference (LHS - RHS)       : {plasma_diff.values[-1]:11.4e}")
        
        print("\n--- Neutral Conservation ---")
        print("    Eq: ddt(N_n) = - Div(Flow_n) - S_iz + S_rec + Puff - Pump + Recycle")
        print(f"    - Div(Flow_n)                : {-total_neutral_flow.values[-1]:11.4e}")
        print(f"    - S_iz                       : {-iz_rate.values[-1]:11.4e}")
        print(f"    + S_rec                      : {recomb_rate.values[-1]:11.4e}")
        print(f"    + Puff                       : {puff_rate.values[-1]:11.4e}")
        print(f"    - Pump                       : {-pump_rate.values[-1]:11.4e}")
        print(f"    + Recycle                    : {target_recycle.values[-1]:11.4e}")
        print(f"    LHS [d/dt(N_n)]              : {ddt_neutral.values[-1]:11.4e}")
        print(f"    RHS [Flows & Sources]        : {neutral_RHS.values[-1]:11.4e}")
        print(f"    Difference (LHS - RHS)       : {neutral_diff.values[-1]:11.4e}")

        print("\n--- Total Conservation ---")
        print("    Eq: ddt(N_tot) = - Div(Flow_tot) + Puff - Pump + Recycle")
        print(f"    - Div(Flow_tot)              : {-total_flow.values[-1]:11.4e}")
        print(f"    + Puff                       : {puff_rate.values[-1]:11.4e}")
        print(f"    - Pump                       : {-pump_rate.values[-1]:11.4e}")
        print(f"    + Recycle                    : {target_recycle.values[-1]:11.4e}")
        print(f"    LHS [d/dt(N_tot)]            : {ddt_total.values[-1]:11.4e}")
        print(f"    RHS [Flows & Sources]        : {total_RHS.values[-1]:11.4e}")
        print(f"    Difference (LHS - RHS)       : {difference.values[-1]:11.4e}")
        
        print(f"\n    True Net Balance             : {true_net_balance.values[-1]:11.4e}")
        print(f"{'='*60}\n")
        
        ax.plot(ds["t"], true_net_balance, label="True Net Balance (Puff-Pump-Loss)", linestyle=":", lw=3, color="magenta", marker = "*")

        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Particle rate [$s^{-1}$]")
        # ax.set_yscale("symlog", linthresh=1e18) 
        # ax.set_ylim([1e19, 1.1e21])
        ax.set_title(f"Particle Balance: {case_name}")
        ax.legend()
        # ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.5)
        
        fig.tight_layout()
        fig.savefig(f"{figures_png_path}/particle_balance_{case_name}.png")
        plt.close(fig)

def run_single_plots():

    parser = build_base_parser()
    args = parser.parse_args()
    case = read_files(args.input)

    ## create output directory

    figures_png_path = args.output + "_figures_png"

    if not os.path.exists(f"./{figures_png_path}"):
        os.makedirs(figures_png_path)

    ## Run
    plot_single_profiles(case, args.region_rad, args.region_pol, args.sepadd, args, figures_png_path)
    # make_plot_diff_coeff(case)
    # plot_Lz_function(case)
    # plot_particle_balance(case, figures_png_path)




from importlib import machinery
from collections import abc
import numpy as np
import matplotlib.pyplot as plt
import shutil
import zipfile
import hashlib
from pathlib import Path
import urllib.request
import glob
import csv
import xhermes
from xhermes.selectors import selector_poloidal, selector_radial, selector_2d
from xhermes.plotting import plot_selection

# ---------------------------------------------------------------------------#
# Settings (change these by hand)
# ---------------------------------------------------------------------------#
verbose = False
plot = True                # save the reconstruction-vs-simulation plots
run_verification = True    # run the reconstruction asserts (set False for table only)
debug = False
rtol = 1e-6
atol = 1e-8
AA_neutral = 2

# Which time step to analyse for EVERY run.
#   -1 = last time step (default)
#   any integer = force every run onto the SAME step index for a fair comparison
timestep = 30

# Core ionising options (src/neutral_boundary.cxx defaults, see BOUT.inp [d] section)
core_ionise_multiplier = 1.0

this_dir = Path(__file__).parent

# ---------------------------------------------------------------------------#
# Multiple input files: every run matching 260617_5e20_corebc*
# ---------------------------------------------------------------------------#
data_dirs = sorted(glob.glob(str(this_dir / "../260617_5e20_corebc*")))
grid_file = f"{this_dir}/../gridfile/CDN_46895_lowresol_260519_nowallpump_5e20.nc"

print("Found runs:")
for d in data_dirs:
    print(f"  {Path(d).name}")
print("\n")


#########################
# Load data
#########################

def load_case(data_dir, unnormalise, timestep=-1):
    """
    Open one Hermes-3 run and return the dataset at the chosen time step.

    unnormalise = False -> normalised units  (used by the reconstruction test)
    unnormalise = True  -> physical SI units (used by the energy-budget table)
    """
    ds = xhermes.open_hermesdataset(
        datapath=f"{data_dir}/BOUT.dmp.*.nc",
        inputfilepath=f"{data_dir}/BOUT.inp",
        gridfilepath=grid_file,
        keep_yboundaries=True,
        keep_xboundaries=True,
        geometry="toroidal",
        unnormalise=unnormalise,
    )

    # Pick the requested time step.
    # Safety guard: a single-time-step dataset may keep t as a length-1 dim,
    # so we select again if it is still present.
    if "t" in ds.sizes:
        ds = ds.isel(t=timestep)

    return ds


def at_boundary(f, i, g):
    """
    Linear interpolation to cell face between domain cell (i) and guard cell (g).
    f = field, i_slice = domain cell index, g_slice = guard cell index.
    """
    return 0.5 * (f[i].values + f[g].values)


def reconstruct_core_ionising(ds):
    """
    Mirrors the core ionising code in neutral_boundary.cxx (lines ~475-589).

    Neutral particle/momentum/energy flowing into the core boundary is removed
    from the neutral population and turned into an ion particle source.
    Momentum/energy only return to the ions if ionising_core_return_mom_energy
    is set (default False). A fraction of the ionisation energy
    (ionising_core_iz_energy_loss, default 13.6 eV) is removed from the
    electron energy as the ionisation cost.

    Works on a NORMALISED dataset (unnormalise=False).
    """

    # Read normalisation and the per-run options straight from this dataset
    Tnorm = ds.metadata["Tnorm"]
    options = ds.options
    val = options["d"]["ionising_core_return_mom_energy"]
    ionising_core_return_mom_energy = str(val).strip().lower() in ("true", "1", "yes")
    ionising_core_iz_energy_loss = options["d"]["ionising_core_iz_energy_loss"]

    if verbose:
        print(f"return_energy_mom (bool) = {ionising_core_return_mom_energy}")
        print(f"iz_energy_loss = {ionising_core_iz_energy_loss}")

    # Select region
    # low: towards core
    x_domain = selector_radial(ds, "boundary_xlow")
    x_guard = selector_radial(ds, "boundary_guard_xlow")
    y_domain = selector_poloidal(ds, "core")
    if verbose:
        print(f"x (boundary low) = {x_domain}, x_guard (boundary guard low) = {x_guard}")
        print(f"y (core) = {y_domain}")

    idx_domain = (x_domain, y_domain)
    idx_guard = (x_guard, y_domain)
    # If you're using neumann boundary condition, values at domain and guard are the same.

    # Select region for diagnostic
    # diagnostic only exist in domain cells, so only needs boundary region
    ds_core = ds.hermes.select_region("boundary", "core")

    Nd_boundary = at_boundary(ds["Nd"], idx_domain, idx_guard)
    Td_boundary = at_boundary(ds["Td"], idx_domain, idx_guard)
    NVd_boundary = at_boundary(ds["NVd"], idx_domain, idx_guard)

    ### Calculate expected values at core boundary

    # cell volume and area at boundary
    dy_boundary = at_boundary(ds["dy"], idx_domain, idx_guard)
    dz_boundary = at_boundary(ds["dz"], idx_domain, idx_guard)
    dpol_boundary = (dy_boundary) / (0.5 * (np.sqrt(ds["g22"][idx_domain].values) + np.sqrt(ds["g22"][idx_guard].values)))
    dtor_boundary = (dz_boundary * 0.5 * (np.sqrt(ds["g_33"][idx_domain].values) + np.sqrt(ds["g_33"][idx_guard].values)))
    da_boundary = dtor_boundary * dpol_boundary
    dv = ds["J"][idx_domain] * ds["dx"][idx_domain] * ds["dy"][idx_domain] * ds["dz"][idx_domain]

    # Mean thermal speed (full one-sided Maxwellian, Stangeby eq. 2.21)
    v_th = np.sqrt(8 * Td_boundary / (np.pi * AA_neutral))

    # Neutral flows reaching the core boundary (Stangeby eq. 2.24: Gamma = n*v_th/4)
    neutral_pflow_to_core = 0.25 * v_th * Nd_boundary * da_boundary
    neutral_mflow_to_core = 0.25 * NVd_boundary * v_th * da_boundary
    neutral_eflow_to_core = 2 * Td_boundary * neutral_pflow_to_core

    # Fraction of neutral flow actually ionised
    ionised_pflow = neutral_pflow_to_core * core_ionise_multiplier
    ionised_mflow = neutral_mflow_to_core * core_ionise_multiplier
    ionised_eflow = neutral_eflow_to_core * core_ionise_multiplier

    # Neutral sinks (always removed regardless of ionising_core_return_mom_energy)
    core_neutral_density_sink = -ionised_pflow / dv
    core_neutral_momentum_sink = -ionised_mflow / dv
    core_neutral_energy_sink = -ionised_eflow / dv

    # Ion particle source always gained
    core_ion_density_source = ionised_pflow / dv

    # Ion momentum/energy source only if removed neutral mom/energy returns as ion mom/energy
    if ionising_core_return_mom_energy:
        core_ion_momentum_source = ionised_mflow / dv
        core_ion_energy_source = ionised_eflow / dv
    else:
        core_ion_momentum_source = np.zeros_like(ionised_mflow)
        core_ion_energy_source = np.zeros_like(ionised_eflow)

    # Electron energy sink: ionisation energy cost paid by the electrons
    if ionising_core_iz_energy_loss != 0.0:
        core_electron_energy_source = -(ionising_core_iz_energy_loss / Tnorm) * ionised_pflow / dv
    else:
        core_electron_energy_source = np.zeros_like(ionised_pflow)

    if verbose:
        print(f"Nd at core = {Nd_boundary}")
        print(f"Td at core = {Td_boundary}")
        print("neutral_pflow_to_core", neutral_pflow_to_core)

    return (
        ds_core,
        core_ion_density_source,
        core_ion_momentum_source,
        core_ion_energy_source,
        core_neutral_density_sink,
        core_neutral_momentum_sink,
        core_neutral_energy_sink,
        core_electron_energy_source,
    )


def plot_result(ax, xvalues, sim_result, calc_result, title):

    ax.plot(xvalues, sim_result, label="simulation result", lw=1, ms=7, marker="o", c="k", )
    ax.plot(xvalues, calc_result, label="calc result", lw=1, ms=7, marker="o", c="r")
    ax.set_ylabel("Value", fontsize=14)
    ax2 = ax.twinx()
    ax2.plot(xvalues, calc_result / sim_result, lw=1, c="blue")
    ax.set_title(title, fontsize=14)
    ax.legend(fontsize=14)


# -------------------------------------------------------------------------------#
### Verification test: reconstruction vs simulation, run for EACH file
# -------------------------------------------------------------------------------#

def verify_core_ionising(data_dir, timestep=-1):
    """Reconstruct the core-ionising sources for one run and check they match
    the simulation diagnostics. Saves a comparison plot per run."""

    name = Path(data_dir).name
    print(f"\n--- Verifying reconstruction: {name} ---")

    # Verification works in normalised units
    ds = load_case(data_dir, unnormalise=False, timestep=timestep)

    (
        ds_core,
        core_ion_density_source,
        core_ion_momentum_source,
        core_ion_energy_source,
        core_neutral_density_sink,
        core_neutral_momentum_sink,
        core_neutral_energy_sink,
        core_electron_energy_source,
    ) = reconstruct_core_ionising(ds)

    if plot:
        fig, ax = plt.subplots(3, 3, figsize=(18, 18))
        plot_result(ax[0, 0], ds_core["theta"].values,
                    ds_core["Nd_core_sink"].values, core_neutral_density_sink,
                    "neutral particle sink at core [m^-3 s^-1]")
        plot_result(ax[0, 1], ds_core["theta"].values,
                    ds_core["NVd_core_sink"].values, core_neutral_momentum_sink,
                    "neutral momentum sink at core [kg m^-2 s^-2]")
        plot_result(ax[0, 2], ds_core["theta"].values,
                    ds_core["Ed_core_sink"].values, core_neutral_energy_sink,
                    "neutral energy sink at core [W m^-3]")
        plot_result(ax[1, 0], ds_core["theta"].values,
                    ds_core["Nd+_core_source"].values, core_ion_density_source,
                    "ion particle source at core [m^-3 s^-1]")
        plot_result(ax[1, 1], ds_core["theta"].values,
                    ds_core["NVd+_core_source"].values, core_ion_momentum_source,
                    "ion momentum source at core [kg m^-2 s^-2]")
        plot_result(ax[1, 2], ds_core["theta"].values,
                    ds_core["Ed+_core_source"].values, core_ion_energy_source,
                    "ion energy source at core [W m^-3]")
        plot_result(ax[2, 0], ds_core["theta"].values,
                    ds_core["Ee_core_source"].values, core_electron_energy_source,
                    "electron energy sink (ionisation loss) at core [W m^-3]")
        fig.tight_layout()
        fig.savefig(f"core_ionising_test_{name}.png")
        plt.close(fig)

    ### Asserts
    np.testing.assert_allclose(core_neutral_density_sink, ds_core["Nd_core_sink"].values,
                               rtol=rtol, atol=atol, err_msg="Neutral particle sink mismatch")
    print("  Neutral particle sink match")
    np.testing.assert_allclose(core_neutral_momentum_sink, ds_core["NVd_core_sink"].values,
                               rtol=rtol, atol=atol, err_msg="Neutral momentum sink mismatch")
    print("  Neutral momentum sink match")
    np.testing.assert_allclose(core_neutral_energy_sink, ds_core["Ed_core_sink"].values,
                               rtol=rtol, atol=atol, err_msg="Neutral energy sink mismatch")
    print("  Neutral energy sink match")
    np.testing.assert_allclose(core_ion_density_source, ds_core["Nd+_core_source"].values,
                               rtol=rtol, atol=atol, err_msg="Ion density source mismatch")
    print("  Ion density source match")
    np.testing.assert_allclose(core_ion_momentum_source, ds_core["NVd+_core_source"].values,
                               rtol=rtol, atol=atol, err_msg="Ion momentum source mismatch")
    print("  Ion momentum source match")
    np.testing.assert_allclose(core_ion_energy_source, ds_core["Ed+_core_source"].values,
                               rtol=rtol, atol=atol, err_msg="Ion energy source mismatch")
    print("  Ion energy source match")
    np.testing.assert_allclose(core_electron_energy_source, ds_core["Ee_core_source"].values,
                               rtol=rtol, atol=atol, err_msg="Electron energy source (ionisation loss) mismatch")
    print("  Electron energy source match")


# -------------------------------------------------------------------------------#
### Energy budget table: compare ionising_core energy flows across all runs
# -------------------------------------------------------------------------------#

def analyse_core_energy_budget(timestep=-1, csv_path="core_energy_budget.csv"):
    """
    Integrate the core-ionising energy diagnostics over the core boundary ring
    for every run, print a comparison table (in kW) and save it as a CSV file.

    Columns (all volume-integrated over the core boundary cells):
      Ion energy return : Ed+_core_source  (energy handed back to the ions)
      e- izloss         : Ee_core_source   (electron ionisation cost, the izloss flag)
      Neutral sink      : Ed_core_sink     (neutral energy removed, always on)
      Net offset        : sum of the three -> net energy added/removed from the system
      e- izloss / sink  : |izloss| / |neutral sink| as a percentage
    """

    print("=" * 110)
    print(f" CORE IONISING ENERGY BUDGET   (time step index = {timestep})")
    print("=" * 110)

    # Markdown table header (renders nicely if you paste it into notes)
    print("| Run | Ion energy return | e⁻ izloss | Neutral sink (always on) | "
          "Net system energy offset | e⁻ izloss / sink |")
    print("| --- | --- | --- | --- | --- | --- |")

    # Collect every row so we can write the CSV at the end
    rows = []

    for data_dir in data_dirs:
        # Short, readable label for the table
        name = Path(data_dir).name.replace("260617_5e20_corebc_", "")

        # Energy budget needs physical units (kW), so unnormalise=True
        ds = load_case(data_dir, unnormalise=True, timestep=timestep)
        ds_core = ds.hermes.select_region("boundary", "core")

        # Cell volume of the core boundary cells [m^3]
        dv_core = ds_core["J"] * ds_core["dx"] * ds_core["dy"] * ds_core["dz"]

        def total_kW(varname):
            """Volume-integrate a [W m^-3] diagnostic over the core ring -> kW."""
            if varname not in ds_core:
                return 0.0
            return float((ds_core[varname] * dv_core).sum().values) / 1e3

        ion_return = total_kW("Ed+_core_source")   # energy returned to ions
        e_izloss = total_kW("Ee_core_source")      # electron ionisation cost (izloss)
        n_sink = total_kW("Ed_core_sink")          # neutral energy removed (always on)
        net = ion_return + e_izloss + n_sink       # net energy offset of the system

        pct = abs(e_izloss) / abs(n_sink) * 100 if n_sink != 0 else 0.0

        print(f"| {name} "
              f"| {ion_return:+.2f} kW "
              f"| {e_izloss:+.2f} kW "
              f"| {n_sink:+.2f} kW "
              f"| **{net:+.2f} kW** "
              f"| {pct:.1f}% |")

        rows.append({
            "run": name,
            "ion_energy_return [kW]": round(ion_return, 4),
            "e_izloss [kW]": round(e_izloss, 4),
            "neutral_sink [kW]": round(n_sink, 4),
            "net_energy [kW]": round(net, 4),
            "e_izloss / neutral_sink [%]": round(pct, 2),
        })

    print("=" * 110 + "\n")

    # Write the table to a CSV file (plain numbers, no kW/% suffixes so it stays numeric)
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved table to {csv_path}\n")


# -------------------------------------------------------------------------------#
### Run
# -------------------------------------------------------------------------------#

if __name__ == "__main__":

    # verify the reconstruction matches the simulation for each run
    # if run_verification:
    #     for data_dir in data_dirs:
    #         verify_core_ionising(data_dir, timestep=timestep)

    # Energy budget comparison table (the mainV result)
    analyse_core_energy_budget(timestep=timestep)

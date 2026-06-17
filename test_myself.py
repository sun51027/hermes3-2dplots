
from importlib import machinery
from collections import abc
import numpy as np
import matplotlib.pyplot as plt
import shutil
import zipfile
import hashlib
from pathlib import Path
import urllib.request
import xhermes
from xhermes.selectors import selector_poloidal, selector_radial, selector_2d
from xhermes.plotting import plot_selection

verbose = False
plot = True
debug = False
rtol = 1e-6
atol = 1e-8
AA_neutral = 2

# Core ionising options (src/neutral_boundary.cxx defaults, see BOUT.inp [d] section)
core_ionise_multiplier = 1.0
# ionising_core_return_mom_energy = False  # default: removed neutral mom/energy does NOT return to ions
# ionising_core_iz_energy_loss_eV = 13.6   # default: hydrogen ground-state ionisation potential

this_dir = Path(__file__).parent
data_dir = "../260617_5e20_corebc_return_True_izloss_True"
output_path = f"{data_dir}/BOUT.dmp.*.nc"
input_file = f"{data_dir}/BOUT.inp"
grid_file = f"{this_dir}/../gridfile/CDN_46895_lowresol_260519_nowallpump_5e20.nc"

#########################
# Load data
#########################

ds = xhermes.open_hermesdataset(
    datapath=output_path,
    inputfilepath=input_file,
    gridfilepath=grid_file,
    keep_yboundaries=True,
    keep_xboundaries=True,
    geometry="toroidal",
    unnormalise=False,
).isel(t=-1)

# ds = ds.load()

# Safety guard: after .isel(t=-1), the t dimension should be gone.
# But if the dataset has a single time step, t may persist as a length-1 dim.
if "t" in ds.sizes:  
    ds = ds.isel(t=-1)

# Extract metadata
m = ds.metadata
Tnorm = m["Tnorm"]
# Extract variables from input file
options = ds.options
# print(options)
# core_ionise_multiplier = options["d"]["core_ionise_multiplier"]
# print(core_ionise_multiplier)
ionising_core_return_mom_energy = options["d"]["ionising_core_return_mom_energy"]
ionising_core_iz_energy_loss = options["d"]["ionising_core_iz_energy_loss"]

print("======== Read from BOUT.inp ==============")
print(f"return_energy_mom (bool) = {ionising_core_return_mom_energy}")
print(f"iz_energy_loss = {ionising_core_iz_energy_loss}")
print("\n")

val = options["d"]["ionising_core_return_mom_energy"]
ionising_core_return_mom_energy = str(val).strip().lower() in ("true", "1", "yes")

if verbose:
    print(m)
    print(f"Topology: {m['topology']}")
    print(f"MXG={m['MXG']}, MYG={m['MYG']}")
    print(f"nxg={m['nxg']}, nyg={m['nyg']}")

def at_boundary(f, i, g):
    """
    Linear interpolation to cell face between domain cell (i) and guard cell (g).
    f = field, i_slice = domain cell index, g_slice = guard cell index.
    """
    return 0.5 * (f[i].values + f[g].values)

def reconstruct_core_ionising():
    """
    Mirrors the core ionising code in neutral_boundary.cxx (lines ~475-589).

    Neutral particle/momentum/energy flowing into the core boundary is removed
    from the neutral population and turned into an ion particle source.
    Momentum/energy only return to the ions if ionising_core_return_mom_energy
    is set (default False). A fraction of the ionisation energy
    (ionising_core_iz_energy_loss, default 13.6 eV) is removed from the
    electron energy as the ionisation cost.
    """

    # Select region
    # low: towards core
    x_domain = selector_radial(ds,"boundary_xlow")
    x_guard = selector_radial(ds,"boundary_guard_xlow")
    y_domain = selector_poloidal(ds, "core")
    if verbose:
        print(f"x (boundary low) = {x_domain}, x_guard (boundary guard low) = {x_guard}")
        print(f"y (core) = {y_domain}")

    idx_domain = (x_domain, y_domain)
    idx_guard = (x_guard, y_domain)
    # If you're using neumann boundary condition, values at domain and guard are the same.
    # print(ds["Nd"][idx_guard].values)

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
    dpol_boundary =   (dy_boundary) / (0.5 * (np.sqrt(ds["g22"][idx_domain].values)+np.sqrt(ds["g22"][idx_guard].values)))
    dtor_boundary =   (dz_boundary * 0.5 * (np.sqrt(ds["g_33"][idx_domain].values)+np.sqrt(ds["g_33"][idx_guard].values)))
    da_boundary = dtor_boundary * dpol_boundary
    dv = ds["J"][idx_domain] * ds["dx"][idx_domain] * ds["dy"][idx_domain] * ds["dz"][idx_domain]

    # Mean thermal speed (full one-sided Maxwellian, Stangeby eq. 2.21)
    v_th = np.sqrt(8 * Td_boundary/ (np.pi * AA_neutral ))

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
        print("neutral_pflow_to_core",neutral_pflow_to_core)

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

    ax.plot(xvalues, sim_result, label = "simulation result", lw=1, ms=7, marker="o", c="k")
    ax.plot(xvalues, calc_result, label = "calc result", lw=1, ms=7, marker="o", c="r")
    ax.set_ylabel("Value")
    ax2 = ax.twinx()
    ax2.plot(xvalues, calc_result/sim_result, lw=1, c="blue")
    ax.set_title(title, fontsize="small")
    ax.legend(fontsize="x-small")

# -------------------------------------------------------------------------------#
### Test
# -------------------------------------------------------------------------------#

(
    ds_core,
    core_ion_density_source,
    core_ion_momentum_source,
    core_ion_energy_source,
    core_neutral_density_sink,
    core_neutral_momentum_sink,
    core_neutral_energy_sink,
    core_electron_energy_source,
) = reconstruct_core_ionising()

if plot:
    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    plot_result(
        ax[0, 0],
        ds_core["theta"].values,
        sim_result = ds_core["Nd_core_sink"].values,
        calc_result = core_neutral_density_sink,
        title = f"neutral particle sink at core [m^-3 s^-1]"
    )
    plot_result(
        ax[0, 1],
        ds_core["theta"].values,
        sim_result = ds_core["NVd_core_sink"].values,
        calc_result = core_neutral_momentum_sink,
        title = f"neutral momentum sink at core [kg m^-2 s^-2]"
    )
    plot_result(
        ax[0, 2],
        ds_core["theta"].values,
        sim_result = ds_core["Ed_core_sink"].values,
        calc_result = core_neutral_energy_sink,
        title = f"neutral energy sink at core [W m^-3]"
    )
    plot_result(
        ax[1, 0],
        ds_core["theta"].values,
        sim_result = ds_core["Nd+_core_source"].values,
        calc_result = core_ion_density_source,
        title = f"ion particle source at core [m^-3 s^-1]"
    )
    plot_result(
        ax[1, 1],
        ds_core["theta"].values,
        sim_result = ds_core["NVd+_core_source"].values,
        calc_result = core_ion_momentum_source,
        title = f"ion momentum source at core [kg m^-2 s^-2]"
    )
    plot_result(
        ax[1, 2],
        ds_core["theta"].values,
        sim_result = ds_core["Ee_core_source"].values,
        calc_result = core_electron_energy_source,
        title = f"electron energy sink (ionisation loss) at core [W m^-3]"
    )
    fig.tight_layout()
    fig.savefig("core_ionisin_test.png")

### Asserts
np.testing.assert_allclose(
    core_neutral_density_sink,
    ds_core["Nd_core_sink"].values,
    rtol=rtol,
    atol=atol,
    err_msg="Neutral particle sink mismatch"
)
print("Neutral particle sink match")
np.testing.assert_allclose(
    core_neutral_momentum_sink,
    ds_core["NVd_core_sink"].values,
    rtol=rtol,
    atol=atol,
    err_msg="Neutral momentum sink mismatch"
)
print("Neutral momentum sink match")
np.testing.assert_allclose(
    core_neutral_energy_sink,
    ds_core["Ed_core_sink"].values,
    rtol=rtol,
    atol=atol,
    err_msg="Neutral energy sink mismatch"
)
print("Neutral energy sink match")
np.testing.assert_allclose(
    core_ion_density_source,
    ds_core["Nd+_core_source"].values,
    rtol=rtol,
    atol=atol,
    err_msg="Ion density source mismatch"
)
print("Ion density source match")
np.testing.assert_allclose(
    core_ion_momentum_source,
    ds_core["NVd+_core_source"].values,
    rtol=rtol,
    atol=atol,
    err_msg="Ion momentum source mismatch"
)
print("Ion momentum source match")
np.testing.assert_allclose(
    core_ion_energy_source,
    ds_core["Ed+_core_source"].values,
    rtol=rtol,
    atol=atol,
    err_msg="Ion energy source mismatch"
)
print("Ion energy source match")
np.testing.assert_allclose(
    core_electron_energy_source,
    ds_core["Ee_core_source"].values,
    rtol=rtol,
    atol=atol,
    err_msg="Electron energy source (ionisation loss) mismatch"
)
print("Electron energy source match")

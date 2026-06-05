
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

this_dir = Path(__file__).parent
data_dir = "../260528_5e20_neumann_nowallpump_coreiz"
output_path = f"{data_dir}/BOUT.dmp.*.nc"
input_file = f"{data_dir}/BOUT.inp"
grid_file = f"{this_dir}/../gridfile/CDN_46895_lowresol_260519_nowallpump_5e20.nc"

#########################
# Load data
#########################

ds = xhermes.open_hermesdataset(
    datapath=output_path,
    # inputfilepath=input_file,
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
    
    ### Calculate expected values at core boundary
    
    # cell volume and area at boundary
    dy_boundary = at_boundary(ds["dy"], idx_domain, idx_guard) 
    dz_boundary = at_boundary(ds["dz"], idx_domain, idx_guard) 
    dpol_boundary =   (dy_boundary) / (0.5 * (np.sqrt(ds["g22"][idx_domain].values)+np.sqrt(ds["g22"][idx_guard].values)))
    dtor_boundary =   (dz_boundary * 0.5 * (np.sqrt(ds["g_33"][idx_domain].values)+np.sqrt(ds["g_33"][idx_guard].values)))
    da_boundary = dtor_boundary * dpol_boundary
    dv = ds["J"][idx_domain] * ds["dx"][idx_domain] * ds["dy"][idx_domain] * ds["dz"][idx_domain]
    
    # Sound speed
    v_th = np.sqrt(8 * Td_boundary/ (np.pi * AA_neutral ))
    
    # Reconstruct expected values at core boundary
    neutral_pflow_to_core = 0.25 * v_th  * Nd_boundary * da_boundary  
    # neutral_mflow_to_core = 0.25 * AA_neutral * Ndcore * v_th**2 * da_boundary 
    neutral_mflow_to_core = 1/2 * Nd_boundary * Td_boundary * da_boundary  
    neutral_eflow_to_core = 2 * Td_boundary * neutral_pflow_to_core

    core_ion_density_source = neutral_pflow_to_core / dv
    core_ion_momentum_source = neutral_mflow_to_core / dv
    core_ion_energy_source = neutral_eflow_to_core /dv

    
    if verbose:
        print(f"Nd at core = {Nd_boundary}")
        print(f"Td at core = {Td_boundary}")
        # print(f"dacore: {dacore}")
        # print("v_th",v_th)
        # print(f"dv = {dv.values}")
        print("neutral_pflow_to_core",neutral_pflow_to_core)
    
    return ds_core, core_ion_density_source, core_ion_momentum_source, core_ion_energy_source 
    


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

ds_core, core_ion_density_source, core_ion_momentum_source, core_ion_energy_source = reconstruct_core_ionising()

if plot:
    fig, ax = plt.subplots(1, 3, figsize=(18, 6))
    plot_result(
        ax[0], 
        ds_core["theta"].values, 
        sim_result = np.abs(ds_core["Nd_core_sink"]).values,
        calc_result = core_ion_density_source,
        title = f"neutral particle sink = ion particle source at core [m^-3 s^-1]"
    )
    plot_result(
        ax[1], 
        ds_core["theta"].values, 
        sim_result = np.abs(ds_core["NVd_core_sink"]).values,
        calc_result = core_ion_momentum_source,
        title = f"neutral momentum sink = ion momentum source at core [kg m^-2 s^-2]"
    )
    plot_result(
        ax[2], 
        ds_core["theta"].values, 
        sim_result = np.abs(ds_core["Ed_core_sink"]).values ,
        calc_result = core_ion_energy_source,
        title = f"neutral energy sink = ion energy source at core [W m^-3]"
    )
    fig.tight_layout()
    fig.savefig("core_ionisin_test.png")

### Asserts
np.testing.assert_allclose(
    core_ion_density_source,
    np.abs(ds_core["Nd_core_sink"].values),
    rtol=rtol,
    atol=atol,
    err_msg="Particle flow mismatch"
)
print("Particle flow match")
np.testing.assert_allclose(
    core_ion_momentum_source,
    np.abs(ds_core["NVd_core_sink"].values),
    rtol=rtol,
    atol=atol,
    err_msg="Momentum flow mismatch"
)
print("Momentum flow match")
np.testing.assert_allclose(
    core_ion_energy_source,
    np.abs(ds_core["Ed_core_sink"].values ),
    rtol=rtol,
    atol=atol,
    err_msg="Energy flow mismatch"
)
print("Energy flow match")

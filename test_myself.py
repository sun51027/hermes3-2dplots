
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

verbose = True
plot = True
debug = False
rtol = 1e-6
atol = 1e-8

this_dir = Path(__file__).parent
data_dir = "../260521_1e21_nowallpump_core_iz"
output_path = f"{data_dir}/BOUT.dmp.*.nc"
input_file = f"{data_dir}/BOUT.inp"
grid_file = f"{this_dir}/../gridfile/CDN_46895_lowresol_260519_nowallpump_1e21.nc"

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
# print(m)
print(f"Topology: {m['topology']}")
print(f"MXG={m['MXG']}, MYG={m['MYG']}")
print(f"nxg={m['nxg']}, nyg={m['nyg']}")

# Select region
# low: towards core
x_domain = selector_radial(ds,"boundary_xlow")
x_guard = selector_radial(ds,"boundary_guard_xlow")
print(f"x_domain = {x_domain}, x_guard = {x_guard}")
y_domain = selector_poloidal(ds, "core")
print(f"y_domain = {y_domain}")

idx_domain = (x_domain, y_domain)
idx_guard = (x_guard, y_domain)
# Because of neumann boundary condition, values at domain and guard are the same.
# print(ds["Nd"][idx_guard].values)

# Select region for diagnostic 
# diagnostic only exist in domain cells, so only needs boundary region
ds_core = ds.hermes.select_region("boundary", "core")

AA_neutral = 2

def at_boundary(f, i, g):
    
    return 0.5 * (f[i].values + f[g].values)

Ndcore = at_boundary(ds["Nd"], idx_domain, idx_guard)
Tdcore = at_boundary(ds["Td"], idx_domain, idx_guard)

# Calculate expected values at core boundary
dy_core = at_boundary(ds["dy"], idx_domain, idx_guard) 
dz_core = at_boundary(ds["dz"], idx_domain, idx_guard) 
dpol_core =   (dy_core) / (0.5 * (np.sqrt(ds["g22"][idx_domain].values)+np.sqrt(ds["g22"][idx_guard].values)))
dtor_core =   (dz_core * 0.5 * (np.sqrt(ds["g_33"][idx_domain].values)+np.sqrt(ds["g_33"][idx_guard].values)))
dacore = dtor_core * dpol_core
v_th = np.sqrt(8 * Tdcore/ (np.pi * AA_neutral ))
dv = ds["J"][idx_domain] * ds["dx"][idx_domain] * ds["dy"][idx_domain] * ds["dz"][idx_domain]
neutral_pflow_to_core = 0.25 * v_th  * Ndcore * dacore
# neutral_mflow_to_core = 0.25 * AA_neutral * Ndcore * v_th**2 * dacore 
neutral_mflow_to_core = 1/2 * Ndcore * Tdcore * dacore 
neutral_eflow_to_core = 2 * Tdcore * neutral_pflow_to_core

if verbose:
    print(f"Nd at core = {Ndcore}")
    print(f"Td at core = {Tdcore}")
    # print(f"dacore: {dacore}")
    # print("v_th",v_th)
    # print(f"dv = {dv.values}")
    print("neutral_pflow_to_core",neutral_pflow_to_core)
    


fig, ax = plt.subplots(1,3,figsize=(18,6))
ax[0].plot(ds_core["theta"].values ,neutral_pflow_to_core, label = "expected")
ax[0].plot(ds_core["theta"].values ,np.abs(ds_core["Nd_core_sink"].values * dv), label = "real")
# print(f"expected values : {neutral_pflux_to_core}")
# print(f"real values : {ds_core["Nd_core_sink"].values * dv}")
ax[0].legend()
ax[0].set_xlabel("theta")
ax[0].set_ylabel("neutral_pflux_to_core")
ax[0].set_title("neutral_pflux_to_core")

ax[1].plot(ds_core["theta"].values ,neutral_mflow_to_core, label = "expected")
ax[1].plot(ds_core["theta"].values ,np.abs(ds_core["NVd_core_sink"].values * dv), label = "real")
# print(f"expected values : {neutral_pflux_to_core}")
# print(f"real values : {ds_core["Nd_core_sink"].values * dv}")
ax[1].legend()
ax[1].set_xlabel("theta")
ax[1].set_ylabel("neutral_mflow_to_core")
ax[1].set_title("neutral_mflow_to_core")

ax[2].plot(ds_core["theta"].values ,neutral_eflow_to_core, label = "expected")
ax[2].plot(ds_core["theta"].values ,np.abs(ds_core["Ed_core_sink"].values * dv), label = "real")
# print(f"expected values : {neutral_pflux_to_core}")
# print(f"real values : {ds_core["Nd_core_sink"].values * dv}")
ax[2].legend()
ax[2].set_xlabel("theta")
ax[2].set_ylabel("neutral_eflow_to_core")
ax[2].set_title("neutral_eflow_to_core")
fig.savefig("neutral_pflux_to_core.png")
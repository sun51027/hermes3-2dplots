
# from boututils.run_wrapper import shell, launch_safe, getmpirun
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


"""
This test reconstructs ion density, momentum, and energy sources due to
core ionising, as well as neutral density, momentum, and energy sinks.

It mirrors the core ionising code in neutral_boundary.cxx (lines 478-553).
The calculation determines the neutral flux reaching the core boundary
and converts it into ion sources and neutral sinks.

The core ionising physics (after git update to core_ionising branch):
  1. At the core boundary (x = xstart, i.e. first domain cell),
     calculate midpoint values between domain cell (i) and guard cell (ig).
  2. v_th = sqrt(8*T / (pi*AA))  -- full mean speed (Stangeby eq. 2.21)
  3. Particle flux: Gamma = 0.25 * v_th * n  (one-sided Maxwellian flux, eq. 2.24)
  4. Momentum flux: Pi = 0.5 * n * T  (kinetic pressure on a wall)
  5. Energy flux: Q = 2*T * particle_flow  (Stangeby eq. 2.24: Q = 2kT * Gamma)
  6. Multiply by face area da, divide by cell volume for volumetric source.

Key changes from previous version:
  - v_th no longer carries the 0.25; it is now the full mean thermal speed.
  - Particle flow: 0.25 * v_th * n * da  (was v_th_old * n * da, identical value).
  - Momentum flow: 0.5 * n * T * da  (was v_th_old^2 * n * AA * da -- now pressure-based).
  - Energy flow: 2*T * particle_flow = 0.5*T*v_th*n*da  (was 2*T*v_th*n*da, which was 4x too large).
"""

verbose = True
plot = True
debug = False
rtol = 1e-6
atol = 1e-8

### Setup
## Test in scratch
# scratch_dir = "/users/jpm590/scratch"
# data_dir = f"{scratch_dir}/260521_1e21_nowallpump_core_iz"
# output_path = f"{data_dir}/BOUT.dmp.*.nc"

### Local setup
this_dir = Path(__file__).parent
data_dir = "../260528_5e20_neumann_nowallpump_coreiz/"
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

# Note: extract_2d_tokamak_geometry() is called automatically by
# open_hermesdataset when geometry="toroidal" — no need to call again.

# .load() forces all lazy (dask) arrays into memory as numpy arrays.
# Without it, data is read from disk on every access. After .load(),
# everything is in RAM and subsequent operations are much faster.
ds = ds.load()

# Safety guard: after .isel(t=-1), the t dimension should be gone.
# But if the dataset has a single time step, t may persist as a length-1 dim.
if "t" in ds.sizes:
    ds = ds.isel(t=-1)

# Extract metadata
m = ds.metadata
print(f"Topology: {m['topology']}")
print(f"MXG={m['MXG']}, MYG={m['MYG']}")
print(f"nxg={m['nxg']}, nyg={m['nyg']}")

### Extract options
options = ds.options

# Core ionising parameters from BOUT.inp [d] section
core_ionise_multiplier = 1.0
# core_ionise_multiplier = options["d"].get("core_ionise_multiplier", 1.0)
AA_neutral = 2  # Mass number of neutral deuterium

if verbose:
    print(f"core_ionise_multiplier = {core_ionise_multiplier}")
    print(f"AA_neutral = {AA_neutral}")

# Geometry and physics data as numpy arrays
params = [
    "J",
    "g22",
    "g_22",
    "g_33",
    "dx",
    "dy",
    "dz",
    "dv",
    "Nd",
    "Nd+",
    "Td",
    "Td+",
    "Vd+",
    "Pd",
]
d = {}
for param in params:
    d[param] = ds[param].values


### Functions
def floor(value, minimum):
    return np.maximum(value, minimum)


def at_boundary(f, i_slice, g_slice):
    """
    Linear interpolation to cell face between domain cell (i) and guard cell (g).
    f = field, i_slice = domain cell index, g_slice = guard cell index.
    """
    return 0.5 * (f[i_slice] + f[g_slice])


def reconstruct_core_ionise():
    """
    Reconstructs ion density, momentum, and energy sources due to core ionising,
    as well as neutral density, momentum, and energy sinks.

    This reproduces neutral_boundary.cxx lines 479-569.

    In C++:
      - mesh->firstX() && mesh->periodicY(mesh->xstart) selects the core boundary
      - i  = indexAt(Nn, mesh->xstart, iy, iz)      -> first domain cell
      - ig = indexAt(Nn, mesh->xstart - 1, iy, iz)   -> guard cell
      - Loop runs for ALL iy (0 to LocalNy), including Y guards

    In Python/xhermes:
      - x = MXG is the first domain cell (boundary_xlow)
      - x = MXG-1 is the guard cell (boundary_guard_xlow)
      - We select the "core" poloidal region to match where periodicY is true

    Returns
    -------
    tuple of:
        ds_core : xarray.Dataset
            The core boundary dataset (for extracting expected values).
        ion_density_source : np.ndarray
            Ion density source [normalised m^-3 s^-1].
        ion_momentum_source : np.ndarray
            Ion momentum source [normalised kg m^-2 s^-2].
        ion_energy_source : np.ndarray
            Ion energy source [normalised W m^-3].
        neutral_density_sink : np.ndarray
            Neutral density sink (negative) [normalised m^-3 s^-1].
        neutral_momentum_sink : np.ndarray
            Neutral momentum sink (negative) [normalised kg m^-2 s^-2].
        neutral_energy_sink : np.ndarray
            Neutral energy sink (negative) [normalised W m^-3].
    """

    ### Indexing
    # Radial: first domain cell and its guard cell
    x_i = selector_radial(ds, "boundary_xlow")       # = MXG = 2
    x_g = selector_radial(ds, "boundary_guard_xlow")  # = MXG - 1 = 1

    # Poloidal: core region only (where periodicY is true)
    y_core = selector_poloidal(ds, "core")

    if verbose:
        print(f"x_i (boundary_xlow) = {x_i}")
        print(f"x_g (boundary_guard_xlow) = {x_g}")
        print(f"y_core = {y_core}")

    # Create index tuples for numpy array slicing
    # These correspond to (x, theta) in the 2D arrays
    i = (x_i, y_core)   # First domain cell at core poloidal positions
    g = (x_g, y_core)   # Guard cell at core poloidal positions

    # Select the core boundary region dataset for extracting expected values
    ds_core = ds.hermes.select_region("boundary", "core")

    ### Calculate midpoint values at the core boundary wall
    # This matches C++ lines 497-498:
    #   nncore = 0.5 * (Nn[ig] + Nn[i])
    #   tncore = 0.5 * (Tn[ig] + Tn[i])
    nncore = at_boundary(d["Nd"], i, g)
    tncore = at_boundary(d["Td"], i, g)
    # print("nncore",nncore)

    ### Mean thermal speed (full, one-dimensional)
    # Stangeby p.69 eq. 2.21: <v> = sqrt(8*T / (pi*m))
    # C++ line 499: v_th = sqrt(8 * tncore / (PI * AAn))
    # Note: the 0.25 factor (one-sided flux = n*<v>/4) is applied in particle_flow below.
    v_th = np.sqrt(8 * tncore / (np.pi * AA_neutral))
    

    ### Calculate radial cell face area [m^2]
    # C++ lines 508-513:
    #   dpol = 0.5*(dy[i]+dy[ig]) / (0.5*(sqrt(g22[i])+sqrt(g22[ig])))
    #   dtor = 0.5*(dz[i]+dz[ig]) * 0.5*(sqrt(g_33[i])+sqrt(g_33[ig]))
    #   da = dpol * dtor
    dpol = (
        0.5 * (d["dy"][i] + d["dy"][g])
        / (0.5 * (np.sqrt(d["g22"][i]) + np.sqrt(d["g22"][g])))
    )
    dtor = (
        0.5 * (d["dz"][i] + d["dz"][g])
        * 0.5 * (np.sqrt(d["g_33"][i]) + np.sqrt(d["g_33"][g]))
    )
    da = dpol * dtor
    # print("da", da)

    ### Cell volume of the domain cell that receives the source
    # C++ lines 483-484:
    #   volume = J(xstart, iy) * dx(xstart, iy) * dy(xstart, iy) * dz(xstart, iy)
    volume = d["J"][i] * d["dx"][i] * d["dy"][i] * d["dz"][i]
    # print(f"volume = {volume}")

    ### Multiplier
    multiplier = core_ionise_multiplier

    ### Particle flow [particles/s]
    # One-sided Maxwellian flux: Gamma = n * <v> / 4 = 0.25 * v_th * n
    # C++ line 511: neutral_particle_flow_to_core = 0.25 * v_th * nncore * dacore
    # C++ line 513: ionise_particle_flow = neutral_particle_flow * multiplier
    neutral_particle_flow = 0.25 * v_th * nncore * da
    ionise_particle_flow = neutral_particle_flow * multiplier

    ### Momentum flow [kg m/s per s = force]
    # Pressure on a wall from a static Maxwellian: Pi = n*T/2
    # C++ line 517: neutral_momentum_flow_to_core = 0.5 * nncore * tncore * dacore
    # C++ line 518: ionise_momentum_flow = neutral_momentum_flow * multiplier
    # Physics note: this is the momentum flux (pressure) to the wall = nT/2 * area.
    # nT is the normalised pressure (Pn), so 0.5*n*T = P/2.
    neutral_momentum_flow = 0.5 * nncore * tncore * da
    ionise_momentum_flow = neutral_momentum_flow * multiplier

    ### Energy flow [W = J/s]
    # C++ (corrected): neutral_energy_flow_to_core = 2.0 * tncore * neutral_particle_flow_to_core
    #
    # Physics: For a static Maxwellian, the energy flux is Q = 2kT * Gamma
    #   where Gamma = n*<v>/4 is the one-sided particle flux (Stangeby eq. 2.24).
    #   So Q = 2*T * particle_flow  (particle_flow already has the 1/4 and area factors).
    #
    # The previous formula 2*T*v_th*n*da was 4x too large because it used
    # the full v_th without the 1/4 factor from the particle flux.
    neutral_energy_flow = 2.0 * tncore * neutral_particle_flow
    ionise_energy_flow = neutral_energy_flow * multiplier

    ### Convert flows to volumetric sources/sinks
    # C++ lines 543-554: source = flow / volume, sink = -flow / volume
    ion_density_source = ionise_particle_flow / volume
    ion_momentum_source = ionise_momentum_flow / volume
    ion_energy_source = ionise_energy_flow / volume

    neutral_density_sink = -ionise_particle_flow / volume
    neutral_momentum_sink = -ionise_momentum_flow / volume
    neutral_energy_sink = -ionise_energy_flow / volume

    if verbose:
        print(f"\n--- Core ionising calculation ---")
        print(f"nncore range: [{nncore.min():.6e}, {nncore.max():.6e}]")
        print(f"tncore range: [{tncore.min():.6e}, {tncore.max():.6e}]")
        print(f"v_th range: [{v_th.min():.6e}, {v_th.max():.6e}]")
        print(f"da range: [{da.min():.6e}, {da.max():.6e}]")
        print(f"volume range: [{volume.min():.6e}, {volume.max():.6e}]")
        print(f"ion_density_source range: [{ion_density_source.min():.6e}, {ion_density_source.max():.6e}]")
        print(f"ion_energy_source range: [{ion_energy_source.min():.6e}, {ion_energy_source.max():.6e}]")

    return (
        ds_core,
        ion_density_source,
        ion_momentum_source,
        ion_energy_source,
        neutral_density_sink,
        neutral_momentum_sink,
        neutral_energy_sink,
    )


def plot_result(ax, xplot, sim_result, calc_result, title):
    ax.plot(xplot, sim_result, label="Simulation", lw=1, ms=7, marker="o", c="k")
    ax.plot(xplot, calc_result, label="Calculated", lw=1, ms=5, marker="x", c="r")
    ax.set_ylabel("Value")
    ax2 = ax.twinx()
    nonzero = sim_result != 0
    ratio = np.where(nonzero, calc_result / sim_result, 1.0)
    ax2.plot(xplot, ratio, lw=1, c="blue", alpha=0.5)
    ax2.set_ylabel("Calc/sim", color="blue")
    ax.set_title(title, fontsize="small")
    ax.legend(fontsize="x-small")


# -------------------------------------------------------------------------------#
### Test
# -------------------------------------------------------------------------------#

### Calculate data
(
    ds_core,
    ion_density_source,
    ion_momentum_source,
    ion_energy_source,
    neutral_density_sink,
    neutral_momentum_sink,
    neutral_energy_sink,
) = reconstruct_core_ionise()

### Extract expected values from simulation output
# These are the diagnostic fields saved by neutral_boundary.cxx::outputVars
expected_ion_density_source = ds_core["Nd+_core_source"].values
expected_ion_momentum_source = ds_core["NVd+_core_source"].values
expected_ion_energy_source = ds_core["Ed+_core_source"].values

expected_neutral_density_sink = ds_core["Nd_core_sink"].values
expected_neutral_momentum_sink = ds_core["NVd_core_sink"].values
expected_neutral_energy_sink = ds_core["Ed_core_sink"].values

### Plotting
if plot:
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), dpi=100)
    fig.suptitle("Core Ionising: Python Reconstruction vs Simulation", fontsize=14)

    # Use theta coordinate for x-axis
    theta_plot = ds_core["theta"].values

    # Row 0: Ion sources
    plot_result(
        axes[0, 0], theta_plot,
        sim_result=expected_ion_density_source,
        calc_result=ion_density_source,
        title="Ion density source (Nd+_core_source)",
    )
    plot_result(
        axes[0, 1], theta_plot,
        sim_result=expected_ion_momentum_source,
        calc_result=ion_momentum_source,
        title="Ion momentum source (NVd+_core_source)",
    )
    plot_result(
        axes[0, 2], theta_plot,
        sim_result=expected_ion_energy_source,
        calc_result=ion_energy_source,
        title="Ion energy source (Ed+_core_source)",
    )

    # Row 1: Neutral sinks
    plot_result(
        axes[1, 0], theta_plot,
        sim_result=expected_neutral_density_sink,
        calc_result=neutral_density_sink,
        title="Neutral density sink (Nd_core_sink)",
    )
    plot_result(
        axes[1, 1], theta_plot,
        sim_result=expected_neutral_momentum_sink,
        calc_result=neutral_momentum_sink,
        title="Neutral momentum sink (NVd_core_sink)",
    )
    plot_result(
        axes[1, 2], theta_plot,
        sim_result=expected_neutral_energy_sink,
        calc_result=neutral_energy_sink,
        title="Neutral energy sink (Ed_core_sink)",
    )

    fig.tight_layout()
    fig.savefig("core_ionising_test.png", bbox_inches="tight")
    if verbose:
        print("\nPlot saved to core_ionising_test.png")

### Asserts
print("\n--- Running assertions ---")

np.testing.assert_allclose(
    ion_density_source,
    expected_ion_density_source,
    rtol=rtol, atol=atol,
    err_msg="Ion density source mismatch",
)
print("✓ Ion density source matches")

np.testing.assert_allclose(
    ion_momentum_source,
    expected_ion_momentum_source,
    rtol=rtol, atol=atol,
    err_msg="Ion momentum source mismatch",
)
print("✓ Ion momentum source matches")

np.testing.assert_allclose(
    ion_energy_source,
    expected_ion_energy_source,
    rtol=rtol, atol=atol,
    err_msg="Ion energy source mismatch",
)
print("✓ Ion energy source matches")

np.testing.assert_allclose(
    neutral_density_sink,
    expected_neutral_density_sink,
    rtol=rtol, atol=atol,
    err_msg="Neutral density sink mismatch",
)
print("✓ Neutral density sink matches")

np.testing.assert_allclose(
    neutral_momentum_sink,
    expected_neutral_momentum_sink,
    rtol=rtol, atol=atol,
    err_msg="Neutral momentum sink mismatch",
)
print("✓ Neutral momentum sink matches")

np.testing.assert_allclose(
    neutral_energy_sink,
    expected_neutral_energy_sink,
    rtol=rtol, atol=atol,
    err_msg="Neutral energy sink mismatch",
)
print("✓ Neutral energy sink matches")

print("\n=== All core ionising tests passed! ===")

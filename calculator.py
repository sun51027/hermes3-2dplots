# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: venv (3.12.3)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## Quick check any physics parameters for a single file

# %%


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys, pathlib, shlex, subprocess
import xbout
import scipy
import xhermes
from xhermes import *


sys.path.append(r"/users/jpm590/2dspace/post-processing/sdtools/")


from hermes3.utils import *


from hermes3.fluxes import *
from hermes3.case_db import *
from hermes3.load import *
from hermes3.named_selections import *
from hermes3.plotting import *
from hermes3.grid_fields import *
from hermes3.accessors import *
from hermes3.selectors import *

# %load_ext autoreload
# %autoreload 2
print("Done")


# %% [markdown]
#

# %%
db = CaseDB(
    case_dir = r"/users/jpm590/scratch/",
    # grid_dir = r"/users/jpm590/2dspace/hermes-3/build-mc-master"
    grid_dir = r"/users/jpm590/neutralrun/hermes-3/master"

)

toload = [
    # dict(name="source_guess", id="260112-cdn-46895-david-param", unnormalise_geom = True, use_xhermes = True, squash = True),
    # dict(name="david_param", id="260105-cdn-46895-david-param", unnormalise_geom = True, use_xhermes = True, squash = True),
    dict(name="nowallpump_2e21", id="260207-cdn-46895-nowallpump_2e21", unnormalise_geom = True, use_xhermes = True, squash = True),
    dict(name="nowallpump_1e21", id="260207-cdn-46895-nowallpump_1e21", unnormalise_geom = True, use_xhermes = True, squash = True),
    dict(name="nowallpump_1e21_decayBC", id="260207-cdn-46895-nowallpump_decayBC_1e21", unnormalise_geom = True, use_xhermes = True, squash = True),
    dict(name="nowallpump_1e21_limiterOFF", id="260207-cdn-46895-nowallpump_limiterOff_1e21", unnormalise_geom = True, use_xhermes = True, squash = True),


    # dict(name="new_grid", id="251210-cdn-46895", unnormalise_geom = True, use_xhermes = True, squash = True),
    # dict(name="original", id="251007-2D-MASTU", unnormalise_geom = True, use_xhermes = True, squash = True)
    

]
cs = {}
for case in toload:
    cs[case["name"]] = db.load_case_2D(case["id"], use_squash = case["squash"], verbose = True)
    cs[case["name"]].extract_2d_tokamak_geometry()

# %%
### Check meta data
m = cs["original"].ds.metadata
# m
ds = cs["original"].ds
# ds



# %% [markdown]
# ### Calculate the total power source
# By `Pd+_src` + `Pe_src` times dv
#
# Energy = 3/2 * pressure * volume 
#
# `P_src` are in per second. So it's Power = 3/2 * P_src * dv

# %%

######
ds = cs["nowallpump_2e21"].ds.isel(t=-1)
total_press = ds["Pd+_src"].values + ds["Pe_src"].values + ds["Pd_src"].values
power2 = 3/2 * (ds["Pd+_src"].values + ds["Pe_src"].values + ds["Pd_src"].values) * ds["dv"].values
power2 = np.sum(power2)
print(f"Total power source = {power2}")


# %%
#######

ds = cs["nowallpump_1e21"].ds.isel(t=-1)

nparticles = (ds["Sd_src"].values) * ds["dv"].values
nparticles = np.sum(nparticles)
print(f"Total number of particles = {nparticles}")
#######

ds = cs["nowallpump_1e21_decayBC"].ds.isel(t=-1)

nparticles = (ds["Sd_src"].values) * ds["dv"].values
nparticles = np.sum(nparticles)
print(f"Total number of particles = {nparticles}")
#######

ds = cs["nowallpump_1e21_limiterOFF"].ds.isel(t=-1)

nparticles = (ds["Sd_src"].values) * ds["dv"].values
nparticles = np.sum(nparticles)
print(f"Total number of particles = {nparticles}")

#####
ds = cs["nowallpump_2e21"].ds.isel(t=-1)
nparticles = (ds["Sd_src"].values) * ds["dv"].values
nparticles = np.sum(nparticles)
print(f"Total number of particles = {nparticles}")

#####

# %% [markdown]
# ## Calculate volume
#
# Get the volume of the region you want.
#
# You could possibly choose ([reference](https://github.com/mikekryjak/sdtools/blob/4242b6f0a55edf66d0b6e6f706bb7dc8515eb143/hermes3/deprecated/stuff.py#L155)): 
# ```python 
#     core = ds.hermesm.select_region("core_edge")
#     sol = ds.hermesm.select_region("sol_edge")
#     pfr = ds.hermesm.select_region("pfr_edge")
#     domain = ds.hermesm.select_region("all_noguards")
#     domain_volume = domain["dv"].values.sum()
# ```
# When you change a grid, the core volume is changed as well, therefore check "core_edge"

# %%
ds = cs["nowallpump_1e21_limiterOFF"].ds.isel(t=-1)

domain = ds.hermesm.select_region("core_edge")
domain_volume = domain["dv"].values.sum()
print(domain_volume)

#####
ds = cs["nowallpump_1e21"].ds.isel(t=-1)

domain = ds.hermesm.select_region("core_edge")
domain_volume = domain["dv"].values.sum()
print(domain_volume)

# %% [markdown]
# ## Check the position of source and pump
# Polygon plot

# %%

ds = cs["nowallpump_1e21"].ds.isel(t=-1)
fig, ax = plt.subplots(1,3 , figsize=(12,6))
ds["Sd+_src" ].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
ds["Sd_src" ].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )

plt.tight_layout()




# %% [markdown]
# ### Find the first core ring
# The fuel is fuelling in the first core ring, so you should find the correct normalised position for fuelling function in input file.
#
# ```python
# source = H(x) * H({first core ring} - x) * pressure_source
# ```
#
# Below `sep` shows the position of the separatrix. `Srad` is the real distance from the separatrix (-: core, +: Sol). Note that `Srad` is in the cell center, while the real separatrix is on the boundary.
#

# %% [markdown]
# To get the normalised core ring position `x`, take `x` by cumulatively summing up the dx (width of each cell). Then normalised it by dividing max value.
#
# The value you should take is the `x` where `Pe_src` is non-zero (1st core ring)

# %%

ds = cs["david_param"].ds.isel(t=-1)
df = get_1d_radial_data(ds, params = ["dx", "Pd_src","Pe_src","Pd+_src", "Sd_src", "Nd_src"], region = "omp")

df["x"] = df["dx"].cumsum()
# x is in the center of the cell
df["x"] = df["x"]/df["x"].max()
df
# use the first 0.08

# %%

ds = cs["source_guess"].ds.isel(t=-1)
df = get_1d_radial_data(ds, params = ["dx", "Pd_src","Pe_src","Pd+_src", "Sd_src", "Nd_src"], region = "omp")

df["x"] = df["dx"].cumsum()
# x is in the center of the cell
df["x"] = df["x"]/df["x"].max()
df
# use the first 0.08

# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
# ---

# %%

print("hello world")
# # %matplotlib inline

# # %matplotlib qt

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
# ## Read files

# %%
db = CaseDB(
    # case_dir = r"/users/jpm590/2dspace/run",
    case_dir = r"/users/jpm590/scratch/",
    grid_dir = r"/users/jpm590/neutralrun/hermes-3/master"
)

toload = [
    # dict(name="original", id="251007-2D-MASTU", unnormalise_geom = True, use_xhermes = True, squash = True)
    dict(name="MAST-U", id="260212-sruntest_limieteroff", unnormalise_geom = True, use_xhermes = True, squash = True)

]
cs = {}
for case in toload:
    cs[case["name"]] = db.load_case_2D(case["id"], use_squash = case["squash"], verbose = True)
    cs[case["name"]].extract_2d_tokamak_geometry()


# %% [markdown]
# ### temperature profiles

# %%
ds = cs["MAST-U"].ds.isel(t=-1)
fig, ax = plt.subplots(1,3 , figsize=(12,6))

ds["Te"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, vmin=1, vmax=100)
ds["Td"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, vmin=1, vmax=100)
ds["Td+"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True, vmin=1, vmax=100)

plt.tight_layout()

# %% [markdown]
# ### Desnity profiles

# %%
ds = cs["MAST-U"].ds.isel(t=-1)
fig, ax = plt.subplots(1,3 , figsize=(12,6))

ds["Ne"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
ds["Nd"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
ds["Nd+"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True,)

plt.tight_layout()

# %% [markdown]
# ### ddt(Desnity) profiles

# %%
ds = cs["MAST-U"].ds.isel(t=-1)
fig, ax = plt.subplots(1,3 , figsize=(12,6))

ds["ddt(Nd)"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
ds["ddt(Nd+)"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
# %% [markdown]
# ### ddt(Pressure) profiles

# %%
ds = cs["MAST-U"].ds.isel(t=-1)
fig, ax = plt.subplots(1,3 , figsize=(12,6))
ds["ddt(Pd)"].bout.polygon(ax = ax[0], cmap = "Spectral_r", antialias = True, logscale = True, )
ds["ddt(Pd+)"].bout.polygon(ax = ax[1], cmap = "Spectral_r", antialias = True, logscale = True, )
ds["ddt(Pe)"].bout.polygon(ax = ax[2], cmap = "Spectral_r", antialias = True, logscale = True, )

plt.tight_layout()

# %% [markdown]
# ## Time evolution animation

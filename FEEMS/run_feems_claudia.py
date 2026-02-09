#!/usr/bin/env python3 
"""
Script version of the FEEMS Getting Started tutorial:
https://github.com/NovembreLab/feems/blob/main/docsrc/notebooks/getting-started.ipynb
"""

import os
import numpy as np, pickle
import matplotlib 
matplotlib.use("Agg")  # needed for HPC/headless plotting
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, MultiPoint
from sklearn.impute import SimpleImputer
import pickle 

os.environ["CARTOPY_DATA_DIR"] = "/home/qvw641/.conda/envs/feems_e/lib/python3.8/site-packages/cartopy/data/shapefiles/natural_earth/physical"
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from pandas_plink import read_plink
from feems import utils, viz, cross_validation, spatial_graph
from feems import SpatialGraph, Viz

# -------------------------
# USER PATHS (EDIT HERE)
# -------------------------
DATA_DIR = "/projects/psg/people/qvw641/ResQ_wolf/FEEMS/data_exclude/" #all samples after final dataset
#DATA_DIR = "/projects/psg/people/qvw641/ResQ_wolf/FEEMS/data_captive/" #excluding captive individuals
PLINK_PREFIX = os.path.join(DATA_DIR, "modernWolf_pruned")  # no extension (.bed/.bim/.fam)
SAMPLES_TSV = os.path.join(DATA_DIR, "samples.tsv")         # must have columns: sample_id, lon, lat
OUT_PREFIX = os.path.join(DATA_DIR, "wolves_run_all")
# -------------------------


# read_plink returns 3 DataFrames in v1
(bim, fam, bed) = read_plink(PLINK_PREFIX)

print("BIM:", bim.head())
print("FAM:", fam.head())
print("BED shape:", bed.shape)  # (n_snps, n_samples)

# transpose to samples x snps
X = bed.T.astype(float)

# impute missing data if any
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
X_imputed = imp.fit_transform(X)

print("n_samples={}, n_snps={}".format(X_imputed.shape[0], X_imputed.shape[1]))

sample_ids = fam.iid.tolist()

# setup graph
coord = np.loadtxt("{}/modernWolf_pruned.coord".format(DATA_DIR))  # sample coordinates
outer = np.loadtxt("{}/modernWolf_pruned.outer".format(DATA_DIR))  # outer coordinates
grid_path = "{}/grid_100.shp".format(DATA_DIR)  # path to discrete global grid for ASIA
#grid_path = "{}/world_triangle_res7.shp".format(DATA_DIR)  # path to discrete global grid for EUROPE

# make boundary
outer, edges, grid, _ = utils.prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=False, 
                                             buffer=0,
                                             outer=outer)
print("done boundary")
#Setup the SpatialGraph object
sp_graph = SpatialGraph(X_imputed, coord, grid, edges, scale_snps=True)

projection = ccrs.EquidistantConic(central_longitude=13, central_latitude=50)
#projection = ccrs.PlateCarree()  # straight frame

print("done projection")
sp_graph.fit(lamb = 2.0)

# save data

#with open("sp_graph.pkl", "wb") as f:
#    pickle.dump(sp_graph, f)

#print("Saved sp_graph.")


# EUROPE
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
ax.set_extent([-12, 40, 35, 72], crs=ccrs.PlateCarree()) #europe only, maybe ax.set_extent([-12, 60, 35, 72], to include part of russia
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=5, 
        obs_node_size=2.5, sample_pt_color="black", obs_node_alpha=1.0,
        cbar_font_size=10)
v.draw_edges(use_weights=True)

# add ocean
ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="white", edgecolor="none", zorder=200)
ax.coastlines(resolution="50m", linewidth=0.6, zorder=201)

v.draw_obs_nodes(use_ids=False)
#v.draw_edge_colorbar()

fig.savefig(f"{OUT_PREFIX}.surfaceEUROPE_dense_nolegend.png", dpi=300)
fig.savefig(f"{OUT_PREFIX}.surfaceEUROPE_dense_nolegend.pdf", dpi=300, rasterized=True)

# EUROPE with legend
#fig = plt.figure(dpi=300)
#ax = fig.add_subplot(1, 1, 1, projection=projection)
#ax.set_extent([-12, 40, 35, 72], crs=ccrs.PlateCarree()) #europe only, maybe ax.set_extent([-12, 60, 35, 72], to include part of russia
#v = Viz(ax, sp_graph, projection=projection, edge_width=.5,
#        edge_alpha=1, edge_zorder=100, sample_pt_size=5,
#        obs_node_size=2.5, sample_pt_color="black", obs_node_alpha=1.0,
#        cbar_font_size=10)
#v.draw_edges(use_weights=True)

# add ocean
#ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="white", edgecolor="none", zorder=200)
#ax.coastlines(resolution="50m", linewidth=0.6, zorder=201)

#v.draw_obs_nodes(use_ids=False)
#v.draw_edge_colorbar()

#fig.savefig(f"{OUT_PREFIX}.surfaceEUROPE_legend.png", dpi=300)
#fig.savefig(f"{OUT_PREFIX}.surfaceEUROPE_legend.pdf", dpi=300, rasterized=True)



# ASIA
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
ax.set_extent([-15, 165, 20, 80], crs=ccrs.PlateCarree()) #eurasia
v = Viz(ax, sp_graph, projection=projection, edge_width=.5,
        edge_alpha=1, edge_zorder=100, sample_pt_size=5,
        obs_node_size=1.5,
        obs_node_alpha=0.4,
        sample_pt_color="black",
        cbar_font_size=10)
v.draw_edges(use_weights=True)

# add ocean
ax.add_feature(cfeature.OCEAN.with_scale("50m"), facecolor="white", edgecolor="none", zorder=200)
ax.coastlines(resolution="50m", linewidth=0.6, zorder=201)

v.draw_obs_nodes(use_ids=False)
v.draw_edge_colorbar()

fig.savefig(f"{OUT_PREFIX}.surfaceASIA_dense.png", dpi=300)
fig.savefig(f"{OUT_PREFIX}.surfaceASIA_dense.pdf", dpi=300, rasterized=True)

print("finish plotting")


# Install required packages
library(dggridR)
library(sf)

# ---- Build a DGGS with finer resolution than grid_100 ----
# FEEMS: grid_100 --> res = 6
# So pick res = 7 (or 8) for more resolution
dggs <- dgconstruct(
  res        = 7,            # > 6 = finer than grid_100
  precision  = 30,
  projection = "ISEA",
  aperture   = 4,
  topology   = "TRIANGLE"
)

# ---- Generate the full Earth grid and save to shapefile ----
# This uses ONLY official dggridR functions:
# dgearthgrid(dggs, savegrid = "path/to/file.shp")
dgearthgrid(
  dggs,
  savegrid = "Files/world_triangle_res7.shp"   # <-- change PATH
)

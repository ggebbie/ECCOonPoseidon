# 1. read 3d monthly average state variables
# 2. map to sigma-1 vertical coordinate
#    - save to mdsio files
# 3. interpolate to regularpoles grid
#    - save as NetCDF
include("state2sigmaregularpoles.jl")

# 1. take 2dstate and translate it to regularpoles grid.
include("state2d2regularpoles.jl")

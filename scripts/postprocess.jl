#= 1. read 3d monthly average state variables
 2. map native state to sigma-1 vertical coordinate
    - save to mdsio files
 3. interpolate to regularpoles grid
    - save as NetCDF 
 4. put state variables on regularpoles grid.
    - save to NetCDF =# 

include("post_process/state2sigmaregularpoles.jl")

# 1. take state and translate it to regularpoles grid.
include("post_process/state2regularpoles.jl")

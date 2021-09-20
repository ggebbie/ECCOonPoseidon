# vertically map θ, S, p to sigma 1 surfaces.
# then horizontally map on regularpoles grid.
# ggebbie, 20-Sep-2021

# put on regularpoles grid, save as NetCDF
include("mdsio2regularpoles.jl")

# put NetCDF on sigma1
include("netcdffiles2sigma1.jl")

# vertically map θ, S, p to sigma 1 surfaces.
# then horizontally map on regularpoles grid.
# ggebbie, 20-Sep-2021

include("state2sigma1.jl")

# put on regularpoles grid, save as NetCDF
# fails because mdsiofiles wants a meta file that 
# doesn't exist. 
include("mdsiofiles2regularpoles")

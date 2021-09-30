# repeated text that configures regularpoles grid
lat,lon = latlon(γ)

# Set up Cartesian grid for interpolation.
# Time for a structure.
λC,λG,ϕC,ϕG,nx,ny,nyarc,nyantarc,farc,iarc,jarc,warc,fantarc,iantarc,jantarc,wantarc =
    factors4regularpoles(γ)

# get standard levels of MITgcm
z = depthlevels(γ)
nz = length(z)

tstart = 1992 + 1/24
tend = 2018
tecco = range(tstart,step=1/12,stop=2018)
nt = length(tecco)

# reading NetCDF attributes
filelog = rundir(expt)*"available_diagnostics.log"

lonatts = Dict("longname" => "Longitude", "units" => "degrees east")
latatts = Dict("longname" => "Latitude", "units" => "degrees north")
depthatts = Dict("longname" => "Depth", "units" => "m")

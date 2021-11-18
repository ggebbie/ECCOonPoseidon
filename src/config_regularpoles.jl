# repeated text that configures regularpoles grid
lat,lon = latlon(γ)

# Set up Cartesian grid for interpolation.

# get a mask
μ =Γ.hFacC[:,1]
μ[findall(μ.>0.0)].=1.0
μ[findall(μ.==0.0)].=NaN

# lon=[i for i=-179.:2.0:179., j=-89.:2.0:89.]
# lat=[j for i=-179.:2.0:179., j=-89.:2.0:89.]

# (f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))
# λ=(lon=lon,lat=lat,f=f,i=i,j=j,w=w);
#df = DataFrame(f=λ.f[:], i=λ.i[:], j=λ.j[:], w=Float32.(λ.w[:]));
#CSV.write("interp_coeffs.csv", df)

# Time for a structure. Interpolation factors go into last two variables: named tuples.
λC,λG,ϕC,ϕG,nx,ny,nyarc,nyantarc,λarc,λantarc = factors4regularpoles(γ)

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

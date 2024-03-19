# repeated text that configures regularpoles grid

# Set up Cartesian grid for interpolation.
lat,lon = latlon(γ)

# get a mask
μ = wet_mask(Γ)

# μ =Γ.hFacC[:,1]
# μ[findall(μ.>0.0)].=1.0
# μ[findall(μ.==0.0)].=NaN

# interpolation parameters
rp_params = factors4regularpoles(γ)
tecco = times_ecco()

# reading NetCDF attributes
if isdir(rundir(expt))
    filelog = rundir(expt)*"available_diagnostics.log"
end

gridatts = grid_attributes()
# lonatts = Dict("longname" => "Longitude", "units" => "degrees east")
# latatts = Dict("longname" => "Latitude", "units" => "degrees north")
# depthatts = Dict("longname" => "Depth", "units" => "m")

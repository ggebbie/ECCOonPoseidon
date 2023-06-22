#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots
import NaNMath as nm

include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = OHC_helper.wet_pts(Γ)
region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(ocean_mask, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

where_dry = findall(x -> isapprox(x, 0), cell_volumes)

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

#read in the first time step of S and θ
expname = "iter129_bulkformula"
nt = 312

println(expname)
filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

pstdgrid = MeshArray(γ,Float32,50)
for k = 1:50, ff = 1:5
    pstdgrid.f[ff, k] .= -pstdz[k]
end

σ_avg = similar(Float32.(cell_volumes)); fill!(σ_avg, 0.0)
σ = MeshArray(γ,Float32,50);
nt = 12*10
@time for tt in 1:nt
    println(tt)
    fnameS = datafilelist_S[tt]
    fnameθS = datafilelist_θ[tt] #potential temp and abs. salinity

    θS = γ.read(diagpath[expname]*fnameθS,MeshArray(γ,Float32,100))

    θ = θS[:, 1:50]; 
    S = θS[:, 51:end]; 

    p₀ = 2000
    for ijk in eachindex(σ)
        ρref = OHC_helper.densityJMD95.(θ.f[ijk],S.f[ijk], pstdgrid[ijk], p₀) .- 1000. #EOS from MITGCM 
        σ.f[ijk] .= ρref
        σ_avg[ijk] .=  σ_avg[ijk] .+ (σ[ijk] ./ nt)
    end
end

for k = 1:50, ff = 1:5
    isdry = isapprox.(cell_volumes[ff, k], 0.0)
    σ_avg[ff, k][isdry] .= NaN
end

write(datadir(expname * "_AVG_THETA_sigma1.data"), σ_avg)
# volume_average_by_depth(σ_avg, cell_volumes, γ)
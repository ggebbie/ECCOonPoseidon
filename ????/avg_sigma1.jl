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

include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ocean_mask = OHC_helper.wet_pts(Γ)
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
expname = "iter0_bulkformula"
nt = 312

println(expname)
filelist = searchdir(diagpath[expname], "state_2d_set1") # first filter for state_3d_set1
datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

σ_avg = similar(Float32.(cell_volumes)); fill!(σ_avg, 0.0)

pstdgrid = MeshArray(γ,Float32,50)
for k = 1:50, ff = 1:5
    pstdgrid.f[ff, k] .= -pstdz[k]
end

SA = MeshArray(γ,Float32,50); 
Θ = MeshArray(γ,Float32,50);
σ = MeshArray(γ,Float32,50);


@time for tt in 1:nt
    println(tt)
    fnameS = datafilelist_S[tt]
    fnameθS = datafilelist_θ[tt] #potential temp and abs. salinity

    θS = γ.read(diagpath[expname]*fnameθS,MeshArray(γ,Float32,100))

    θ = θS[:, 1:50]; 
    S = θS[:, 51:end]; 

    lon=0.;lat=30.
    p₀ = 1000
    for ijk in eachindex(σ)
        SA.f[ijk] .= gsw_sa_from_sp.(S.f[ijk],pstdgrid.f[ijk],lon,lat)
        Θ.f[ijk] .= gsw_ct_from_pt.(SA.f[ijk],θ.f[ijk])
        σ.f[ijk] .= gsw_rho.(SA.f[ijk],Θ.f[ijk],p₀) .- 1000.
        σ_avg[ijk] .=  σ_avg[ijk] .+ (σ[ijk] ./ nt)
    end
end

for k = 1:50, ff = 1:5
    isdry = isapprox.(cell_volumes[ff, k], 0.0)
    σ_avg[ff, k][isdry] .= NaN
end

write(datadir(expname * "_AVG_THETA_sigma1.data"), σ_avg)
# volume_average_by_depth(σ_avg, cell_volumes, γ)
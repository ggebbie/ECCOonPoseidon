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
@pyimport cmocean.cm as cmo
@pyimport seaborn as sns;
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
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


σ = MeshArray(γ,Float32,50);
θ_list = []; S_list = []; σ_list = []; P_list = []
tt = 1

fnameS = datafilelist_S[tt]
fnameθS = datafilelist_θ[tt] #potential temp and abs. salinity

θS = γ.read(diagpath[expname]*fnameθS,MeshArray(γ,Float32,100))

θ = θS[:, 1:50]; 
S = θS[:, 51:end]; 
p₀ = 2000

for ijk in eachindex(σ)
    ρref = OHC_helper.densityJMD95.(θ.f[ijk],S.f[ijk], pstdgrid[ijk], p₀) .- 1000. #EOS from MITGCM 
    σ.f[ijk] .= ρref
    isdry = isapprox.(cell_volumes[ijk], 0.0)

    θ_list = vcat(θ_list, vec(θ.f[ijk][(!).(isdry)]))
    S_list = vcat(S_list, vec(S.f[ijk][(!).(isdry)]))
    σ_list = vcat(σ_list, vec(σ.f[ijk][(!).(isdry)]))
    P_list = vcat(P_list, vec(pstdgrid.f[ijk][(!).(isdry)]))

end

S_min, S_max = extrema(S_list)
θ_min, θ_max = extrema(θ_list)
S_levels = S_min:0.1:S_max; 

θ_levels = θ_min:0.1:θ_max

S_grid = S_levels' .* ones(length(θ_levels))
θ_grid = ones(length(S_levels))' .* θ_levels

σgrid = OHC_helper.densityJMD95.(θ_grid,S_grid, p₀, p₀) .- 1000. #EOS from MITGCM 

fig, ax = plt.subplots(1,1, figsize = (20, 12.5))
CS = ax.contour(S_grid, θ_grid, σgrid, colors = "black", linewidths = 0.75)
ax.clabel(CS, fontsize=15, inline=true, fmt = "%.1f")
ax.scatter(S_list, θ_list, c = "black", alpha = 0.3)
ax.set_xlabel("Practical Salinity"); ax.set_ylabel("Potential Temperature")
ax.set_title("Pacific Ocean T-S Diagram")
fig
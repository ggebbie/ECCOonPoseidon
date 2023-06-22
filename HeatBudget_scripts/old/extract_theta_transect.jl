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
gr()
import NaNMath as nm
using .OHC_helper
import PyPlot as plt
using PyCall


pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
λ[findall(λ .<= 0) ] = λ[findall(λ .<= 0) ] .+ 360

area = readarea(γ)
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
labels_L = [L"\theta^{129}", L"\theta^{\Delta F}",L"\theta^{\Delta T}", L"\theta^{0}"]
tecco = 1992+1/24:1/12:2018

ocean_mask = wet_pts(Γ)
region = "NPAC30"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

ϕ_mask_Inf = ϕ .* PAC_msk; 
ϕ_mask_Inf[findall(ϕ_mask_Inf .== 0 )] .= Inf
ϕ_mask_min = minimum(ϕ_mask_Inf)

area = readarea(γ)
area_masked = area .* PAC_msk
area_sum = sum(area_masked)

ϕ_S = (ϕ .≈ ϕ_mask_min)
ϕ_S = ϕ_S .* PAC_msk


Zs  = MeshArray(γ,Float32,50); fill!(Zs, 0.0)
for ij in eachindex(Zs)
    Zs.f[ij] .= z[ij[2]]
end

LON  = MeshArray(γ,Float32,50); fill!(LON, 0.0)
for ij in eachindex(Zs)
    LON.f[ij] .= λ[ij[1]]
end
ϕS_big = MeshArray(γ,Float32,50); fill!(ϕS_big, 0.0)
for ij in eachindex(Zs)
    ϕS_big.f[ij] .= ϕ_S[ij[1]]
end
where_θs = findall(ϕS_big[:, lvls] .!= 0 )

function extract_θs(ds, where_θs)
    mask_ds = ds[where_θs]
    vec_list = []
    for ij in eachindex(mask_ds)
        append!(vec_list, vec(mask_ds[ij]))
    end

    return Float32.(vec_list)
end

LON_coords = extract_θs(LON[:, lvls], where_θs)
Z_coords = extract_θs(Zs[:, lvls], where_θs)

expname = "iter129_bulkformula"

filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
bot_lev = lvls[end]; top_lev = lvls[1]    
sθ̄ = []; AdvH = Float32[]; AdvR = Float32[];DiffH = Float32[]; DiffR = Float32[]
tt = 212

fnameθ = datafilelist_θ[tt]
fnameS = datafilelist_S[tt]

sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)

sθ_coords = extract_θs(sθ[:, lvls], where_θs)


fig,axs=plt.subplots(1,1, figsize = (25, 7.5))
ctf = axs.tricontourf(LON_coords, Z_coords, sθ_coords)
cbar = fig.colorbar(ctf)
fig

filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
fname_uvw = datafilelist_uvw[tt]
UVW = γ.read(diagpath[expname]*fname_uvw,MeshArray(γ,Float32,150))
U = UVW[:, 1:50] .* Γ.hFacC; V = UVW[:, 51:100] .* Γ.hFacC;
U = U; V = V
Utr, Vtr = UVtoTransport(U, V, Γ); 
Utr = Utr[:, lvls] 
Vtr = Vtr[:, lvls] 
Etr, Ntr = rotate_uv(Utr, Vtr, Γ);

fig,axs=plt.subplots(1,1, figsize = (25, 7.5))
ctf = axs.tricontourf(LON_coords, Z_coords, 1e-6 .* (extract_θs(Ntr, where_θs)), cmap= "bwr")
cbar = fig.colorbar([ctf][1], extend = "both")
fig



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

@pyimport seaborn as sns;
@pyimport pandas as pd;
colors =  sns.color_palette("deep")[1:4]
sns.set_theme(context = "talk", font_scale = 1.0,
              palette = sns.color_palette("deep"));#sns.set_context("talk")
pygui(false)
cm = pyimport("cmocean.cm");colorway = cm.balance;
include(srcdir("config_exp.jl"))
(ϕ,λ) = latlonC(γ)
area = readarea(γ)
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
labels_L = [L"\theta^{129}", L"\theta^{\Delta F}",L"\theta^{\Delta T}", L"\theta^{0}"]
tecco = 1992+1/24:1/12:2018

ocean_mask = wet_pts(Γ)
region = "NPAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not full")
suffix = "2to3"
uplvl = -2e3; botlvl = -3e3
lvls = findall( botlvl .<= z[:].<= uplvl)

ϕ_mask_Inf = ϕ .* PAC_msk; ϕ_mask_Inf[findall(ϕ_mask_Inf .== 0 )] .= Inf
ϕ_mask_min = minimum(ϕ_mask_Inf)

ϕ_mask_mInf = ϕ .* PAC_msk; ϕ_mask_mInf[findall(ϕ_mask_Inf .== 0 )] .= -Inf
ϕ_mask_max= maximum(ϕ_mask_mInf)


area = readarea(γ)

area_masked = area .* PAC_msk
area_sum = sum(area_masked)

ϕ_S = (ϕ .≈ ϕ_mask_min)
ϕ_S = ϕ_S .* PAC_msk
Area_S = MeshArray(γ,Float32,50)
cell_depths = get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
for ff = 1:5, k = 1:50
    Area_S.f[ff, k] .=  (Γ.DXG.f[ff] .* ϕ_S.f[ff]) .* cell_depths[ff, k] .* (k ∈ lvls)
end
Area_S_sum = sum(Area_S)

ϕ_N = (ϕ .≈ ϕ_mask_max)
ϕ_N = ϕ_N .* PAC_msk


expname = "iter0_bulkformula"
filelist = searchdir(diagpath[expname],"trsp_3d_set2") # first filter for state_3d_set1
datafilelist_H  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"trsp_3d_set3") # first filter for state_3d_set1
datafilelist_R  = filter(x -> occursin("data",x),filelist) # second filter for "data"
nt = length(datafilelist_R)
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"

smush_depths = smush(cell_depths); 
cell_volumes = get_cell_volumes(area, cell_depths);

smush_depths[findall(smush_depths .== 0)] = Inf; 
inv_depths = 1 ./ smush_depths

sθ̄ = []; AdvH = []; AdvR = [];AdvH2 = [];
masked_volume = cell_volumes[:, lvls]
@time for tt=1:10
    fnameH = datafilelist_H[tt]
    fnameR = datafilelist_R[tt]
    fnameθ = datafilelist_θ[tt]
    fnameS = datafilelist_S[tt]

    sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)

    dθλ = γ.read(diagpath[expname]*fnameH,MeshArray(γ,Float32,200))
    uθ = dθλ[:, 151:200]
    vθ = dθλ[:, 151:200]
    Nθ = MeshArray(γ,Float32,50)
    Eθ, Nθ = rotate_uv(uθ, vθ, Γ)
    uvθ = MeshArray(γ,Float32,50)
    wθ_conv = MeshArray(γ,Float32,50)
    fill!(wθ_conv, 0.0)
    calc_UV_conv3D!(uθ, vθ, uvθ);
    wθ = γ.read(diagpath[expname]*fnameR,MeshArray(γ,Float32,50))
    @inbounds wθ_conv.f[ :, 1:49] .= @views wθ.f[:, 1:49] .- wθ.f[:, 2:50]
    push!(AdvH2, sum(Nθ[:, lvls] .* ϕ_S) / sum(masked_volume))

    push!(AdvH, sum(uvθ[:, lvls] .* PAC_msk) / sum(masked_volume))
    push!(AdvR, sum(wθ_conv[:, lvls] .* PAC_msk)  / sum(masked_volume))
    push!(sθ̄, volume_mean(sθ[:, lvls]; weights=masked_volume))
end

dθdt =  Float32.((-AdvH - AdvR))

θ_approx =  cumsum(cat(sθ̄[1], -dθdt .* 2.628e+6, dims = 1))
plot(tecco[1:10], θ_approx[1:end-1], label = "Advective fluxes approximation", 
ylabel = "Temperature", xlabel = "Time")
plot!(tecco[1:10], sθ̄[1:10], label = "true θ")



weighted_means = randn((11, 15))
get_range(x) = maximum(x) - minimum(x)
for i in 1:size(weighted_means)[1]    
    println(weighted_means[i, :])
end

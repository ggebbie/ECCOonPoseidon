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

area = readarea(γ)
area_masked = area .* PAC_msk
area_sum = sum(area_masked)

ϕ_S = (ϕ .≈ ϕ_mask_min)
ϕ_S = ϕ_S .* PAC_msk

cell_depths = get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
smush_depths = smush(cell_depths); 
cell_volumes = get_cell_volumes(area, cell_depths);

smush_depths[findall(smush_depths .== 0)] = Inf; 
inv_depths = 1 ./ smush_depths

Advective_Heat_Budgets = Dict()
masked_volume = cell_volumes[:, lvls]

@time for expname in keys(shortnames)
    Advective_Heat_Budgets[expname] = Dict()
    filelist = searchdir(diagpath[expname],"trsp_3d_set2") # first filter for state_3d_set1
    datafilelist_H  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"trsp_3d_set3") # first filter for state_3d_set1
    datafilelist_R  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist_R)
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    bot_lev = lvls[end]; top_lev = lvls[1]    
    sθ̄ = []; AdvH = []; AdvR = [];AdvH2 = [];
    @time for tt=1:nt
        fnameH = datafilelist_H[tt]
        fnameR = datafilelist_R[tt]
        fnameθ = datafilelist_θ[tt]
        fnameS = datafilelist_S[tt]

        sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)

        dθλ = γ.read(diagpath[expname]*fnameH,MeshArray(γ,Float32,200))
        uθ = dθλ[:, 101:150]; vθ = dθλ[:, 151:200]
        wθ = γ.read(diagpath[expname]*fnameR,MeshArray(γ,Float32,50))

        Eθ, Nθ = rotate_uv(uθ, vθ, Γ)
        wθ_conv = wθ[:, bot_lev+1] .- wθ[:, top_lev] 

        push!(AdvH, sum(Nθ[:, lvls] .* ϕ_S) / sum(masked_volume))
        push!(AdvR, sum(wθ_conv .* PAC_msk)  / sum(masked_volume))
        push!(sθ̄, volume_mean(sθ[:, lvls]; weights=masked_volume))
    end
    AdvH_sum = cumsum(cat(0.0, AdvH .* 2.628e+6, dims = 1))
    AdvR_sum = cumsum(cat(0.0, AdvR .* 2.628e+6, dims = 1))
    Advective_Heat_Budgets[expname]["AdvH"] = AdvH_sum
    Advective_Heat_Budgets[expname]["AdvR"] = AdvR_sum
    Advective_Heat_Budgets[expname]["θ"] = sθ̄
end
AdvH = Advective_Heat_Budgets["iter129_bulkformula"]["AdvH"]
AdvR = Advective_Heat_Budgets["iter129_bulkformula"]["AdvR"]
θ = Advective_Heat_Budgets[expname]["θ"]
jldsave(datadir("HeatBudgetTrue_NPAC_2to3.jld2"); Advective_Heat_Budgets)

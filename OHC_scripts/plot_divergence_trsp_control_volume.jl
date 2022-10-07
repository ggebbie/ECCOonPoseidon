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
#this should actually correspond to the latitude just before 

ϕ_S = (ϕ .≈ ϕ_mask_min)
ϕ_S = ϕ_S .* PAC_msk
Area_S = MeshArray(γ,Float32,50); fill!(Area_S, 0.0)
cell_depths = get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
for ff = 1:5, k = 1:50
    Area_S.f[ff, k] .=  (Γ.DXG.f[ff] .* ϕ_S.f[ff]) .* Γ.DRF[k] 
end
Area_S_sum = sum(Area_S)

area = readarea(γ)
area_masked = area .* PAC_msk; area_sum = sum(area_masked)

smush_depths = smush(cell_depths); 
cell_volumes = get_cell_volumes(area, cell_depths);

smush_depths[findall(smush_depths .== 0)] = Inf; 
inv_depths = 1 ./ smush_depths

masked_volume = cell_volumes[:, lvls]
vol_masked_top = cell_volumes[:, lvls[1]]
vol_masked_bot = cell_volumes[:, lvls[end]]
vol_masked_south_masked = cell_volumes[:, lvls] .* ϕ_S

Advective_Heat_Budgets = Dict()

@time for expname in keys(shortnames)
    Advective_Heat_Budgets[expname] = Dict()
    filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
    datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    nt = length(datafilelist_uvw)
    sθ̄ = []
    bot_lev = lvls[end]; top_lev = lvls[1]
    V_South = []; W_Top = []; W_Bot = []; 
    θ_South = []; θ_Top = []; θ_Bot = []; 
    @time for tt=1:nt
        fname_uvw = datafilelist_uvw[tt]
        fnameS = datafilelist_S[tt]
        fnameθ = datafilelist_θ[tt]
        UVW = γ.read(diagpath[expname]*fname_uvw,MeshArray(γ,Float32,150))
        sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)

        u = UVW[:, 1:50]; u = u .* PAC_msk
        v = UVW[:, 51:100]; v = v .* PAC_msk
        w_top = UVW[:, 100 + top_lev]; 
        w_bot = UVW[:, 100 + bot_lev + 1]; 

        Utr, Vtr = UVtoTransport(u, v, Γ)
        Etr, Ntr = rotate_uv(Utr, Vtr, Γ)

        θ_bot = sum(sθ[:, bot_lev] .* vol_masked_bot)/sum(vol_masked_bot)
        θ_top = sum(sθ[:, top_lev] .* vol_masked_top)/sum(vol_masked_top)

        w_bot = sum(w_bot .* area_masked); w_top = sum(w_top .* area_masked)
        v_south = sum( Ntr[:, lvls] .* ϕ_S);
        θ_south = sum( sθ[:, lvls] .* vol_masked_south_masked) / sum(vol_masked_south_masked);
        push!(V_South, v_south); push!(W_Top, w_top); push!(W_Bot, w_bot); 
        push!(θ_South, θ_south); push!(θ_Top, θ_top); push!(θ_Bot, θ_bot); 
        push!(sθ̄, volume_mean(sθ[:, lvls]; weights = masked_volume))
    end

    AdvH = (V_South.*θ_South) ./  sum(masked_volume)
    AdvR = (W_Bot.*θ_Bot .- W_Top.*θ_Top) ./  sum(masked_volume)

    AdvH_sum = cumsum(cat(0.0, AdvH .* 2.628e+6, dims = 1))
    AdvR_sum = cumsum(cat(0.0, AdvR .* 2.628e+6, dims = 1))
    Advective_Heat_Budgets[expname]["AdvH"] = AdvH_sum
    Advective_Heat_Budgets[expname]["AdvR"] = AdvR_sum
    Advective_Heat_Budgets[expname]["θ"] = sθ̄
end
jldsave(datadir("HeatBudgetApprox_NPAC_2to3.jld2"); Advective_Heat_Budgets)
AdvH = Advective_Heat_Budgets["iter129_bulkformula"]["AdvH"]
AdvR = Advective_Heat_Budgets["iter129_bulkformula"]["AdvR"]
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
    Area_S.f[ff, k] .=  (Γ.DXC.f[ff] .* ϕ_S.f[ff]) .* cell_depths[ff, k] .* (k ∈ lvls)
end
Area_S_sum = sum(Area_S)

expname = "iter129_bulkformula"
filelist = searchdir(diagpath[expname],"trsp_3d_set1") # first filter for state_3d_set1
datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

smush_depths = smush(cell_depths); 
cell_volumes = get_cell_volumes(area, cell_depths);

smush_depths[findall(smush_depths .== 0)] = Inf; 
inv_depths = 1 ./ smush_depths


nt = length(datafilelist_uvw)
Vdθdy = []; Wdθdz_top = []; Wdθdz_bot = []
sθ̄ = []
bot_lev = lvls[end]; top_lev = lvls[1]
masked_volume = cell_volumes[:, lvls]
V_South = []; W_Top = []; W_Bot = []; 

θ_South = []; θ_Top = []; θ_Bot = []; 
θy_South =[];θz_Top = []; θz_Bot = []; 
area_depth_mask = Γ.hFacC.* PAC_msk;
area_depth_mask[findall(area_depth_mask .> 0)] = 1
area_depth_mask[findall(area_depth_mask .<= 0)] = 0

area_masked_top = area_masked .* area_depth_mask[:, lvls[1]]
area_masked_bot = area_masked .* area_depth_mask[:, lvls[end]]

vol_masked_top = cell_volumes[:, lvls[1]]
vol_masked_bot = cell_volumes[:, lvls[end]]
vol_masked_south = cell_volumes .* ϕ_S

Vθ_South = []; Wθ_Top = []; Wθ_Bot = []; 

@time for tt=1:nt
    fname_uvw = datafilelist_uvw[tt]
    fnameS = datafilelist_S[tt]
    fnameθ = datafilelist_θ[tt]
    UVW = γ.read(diagpath[expname]*fname_uvw,MeshArray(γ,Float32,150))
    sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)
    u = UVW[:, 1:50] 
    v = UVW[:, 51:100] 

    E, N = rotate_uv(u, v, Γ)
    N = N .* PAC_msk
    w_top = (UVW[:, 100 + top_lev] .+ UVW[:, 100 + top_lev-1]) ./2 ; 
    w_bot =(UVW[:, 100 + bot_lev] .+ UVW[:, 100 + bot_lev+1]) ./2 ; 

    sθx = MeshArray(γ,Float32,50); sθy = MeshArray(γ,Float32,50)

    # for k in 1:50
    #     tempx, tempy = gradient(sθ[:, k], Γ)
    #     sθx.f[:, k] .= tempx.f; sθy.f[:, k] .= tempy.f
    # end

    # sθtopdiff = (sθ[:, top_lev-1] .- sθ[:, top_lev]) ./ (z[top_lev-1] - z[top_lev])
    # sθbotdiff = (sθ[:, bot_lev+1] .- sθ[:, bot_lev]) ./ (z[bot_lev + 1] - z[bot_lev])
    # θ_bot_interp = (sθ[:, bot_lev] .+ sθ[:, bot_lev + 1]) ./ 2
    # θ_top_interp = (sθ[:, top_lev - 1] .+ sθ[:,top_lev]) ./ 2
    # wθ_bot = w_bot.*θ_bot_interp
    # wθ_bot = sum(wθ_bot .* area_masked_bot)/sum(area_masked_bot)
    # wθ_top = w_bot.*θ_top_interp
    # wθ_top = sum(wθ_top.* area_masked_top)/sum(area_masked_top)
    # vθ_south = N .* sθ
    # vθ_south = sum(vθ_south .* Area_S)/ Area_S_sum;

    θ_bot = sum(sθ[:, bot_lev] .* vol_masked_bot)/sum(vol_masked_bot)
    θ_top = sum(sθ[:, top_lev] .* vol_masked_top)/sum(vol_masked_top)

    w_bot = sum(w_bot.* area_masked_bot)
    w_top = sum(w_top.* area_masked_top)
    v_south = sum( N .* Area_S);
    θ_south = sum( sθ .* vol_masked_south) / sum(vol_masked_south);
    push!(V_South, v_south); push!(W_Top, w_top); push!(W_Bot, w_bot); 
    push!(θ_South, θ_south); push!(θ_Top, θ_top); push!(θ_Bot, θ_bot); 
    # push!(Vθ_South, vθ_south); push!(Wθ_Top, wθ_top); push!(Wθ_Bot, wθ_bot); 


    push!(sθ̄, volume_mean(sθ[:, lvls]; weights = masked_volume))


    # push!(V_south, sum(v .* Area_S) / Area_S_sum)
end

# dθdt = (V_South.*θ_South) .+ (W_Bot.*θ_Bot) .- (W_Top.*θ_Top)
AdvH = (V_South.*θ_South) ./  (sum(cell_volumes[:, lvls]))
AdvR = (W_Bot.*θ_Bot .- W_Top.*(θ_Top)) ./  (sum(cell_volumes[:, lvls]))
dθdt =  AdvH .+ AdvR

# AdvH2 = (Vθ_South.*Area_S_sum) ./  (sum(cell_volumes[:, lvls]))
# AdvR2 = ((Wθ_Bot.*sum(area_masked_bot)) .- (Wθ_Top.*sum(area_masked_top))) ./  (sum(cell_volumes[:, lvls]))
# dθdt2 =  AdvH2 .+ AdvR2

θ_approx =  cumsum(cat(0.0, dθdt .* 2.628e+6, dims = 1))
plot(tecco, θ_approx[1:end-1] ./ mean(θ_approx[1:end-1]) , label = "Advective fluxes approximation", 
ylabel = "θ / mean(θ)", xlabel = "Time")
plot!(tecco, (sθ̄ .- sθ̄[1]) ./ mean(sθ̄.- sθ̄[1]), label = "true θ")

plot(cumsum(AdvH.* 2.628e+6))
plot(cumsum(AdvR.* 2.628e+6))
dsθ̄dt = (sθ̄[2:end] - sθ̄[1:end-1]) ./ 2.628e+6
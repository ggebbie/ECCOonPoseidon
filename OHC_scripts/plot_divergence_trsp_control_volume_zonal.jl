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
region = "NPAC30"; 
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

ϕ_S = (ϕ .== ϕ_mask_min) 
ϕ_S = ϕ_S .* PAC_msk
Area_S = MeshArray(γ,Float32,50)
for ff = 1:5, k = 1:50
    Area_S.f[ff, k] .=  (Γ.DXG.f[ff] .* ϕ_S.f[ff]) .* Γ.DRF[k] .* (k ∈ lvls)
end
Area_S_sum = sum(Area_S)

expname = "iter129_bulkformula"
fileroot = "trsp_3d_set1"
filelist = searchdir(diagpath[expname],fileroot) # first filter for state_3d_set1
datafilelist_uvw  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
datafilelist_S  = filter(x -> occursin("data",x),filelist) # second filter for "data"
filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

cell_depths = get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = get_cell_volumes(area, cell_depths);

smush_depths = smush(cell_depths); 
smush_depths[findall(smush_depths .== 0)] = Inf; 
inv_depths = 1 ./ smush_depths


nt = length(datafilelist_uvw)
Vdθdy = []; Wdθdz_top = []; Wdθdz_bot = []; θ̄ = []
bot_lev = lvls[end]; top_lev = lvls[1]
masked_volumes = cell_volumes[:, lvls]

X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)
idx30 = findfirst( x-> x > 30, X)
area_mask = area .* PAC_msk

ΔY_mask = OHC_helper.ma_zonal_sum( Γ.DYG .* area_mask) ./ OHC_helper.ma_zonal_sum(area_mask)
ΔY_mask[isnan.(ΔY_mask)] .= 0.0
@time for tt=1:nt
    fname_uvw = datafilelist_uvw[tt]
    fnameS = datafilelist_S[tt]
    fnameθ = datafilelist_θ[tt]

    UVW = γ.read(diagpath[expname]*fname_uvw,MeshArray(γ,Float32,150))
    sθ = extract_sθ(expname,diagpath, γ, fnameS, fnameθ, inv_depths)

    v = UVW[:, 51:100]; 
    w = UVW[:, 101:150];

    v̄ = ma_zonal_avg(v, cell_volumes) 
    w̄  = ma_zonal_avg(w, cell_volumes) 
    s̄θ̄ = ma_zonal_avg(sθ, cell_volumes) 
    s̄θ̄dz = (s̄θ̄[1:end-1, :] .- s̄θ̄[2:end, :]) ./ (z[1:end-1] .- z[2:end])
    #convert from ∘/latitude to ∘/m
    s̄θ̄dy =  inv(110574) .* (s̄θ̄[:, 2:end] .- s̄θ̄[:, 1:end-1]) .* inv.(X[2:end]' .- X[1:end-1]')
    vs̄θ̄dy = s̄θ̄dy[lvls, idx30 + 1] .* mean(v̄[lvls, idx30:idx30+1], dims = 2)

    vs̄θ̄dy_avg = sum( vs̄θ̄dy .* ΔzF[lvls]) / sum(ΔzF[lvls])

    w_interp = (w̄[1:end-1, :] .+ w̄[2:end, :]) ./ 2
    ws̄θ̄dz = s̄θ̄dz .* w_interp; ws̄θ̄dz[isnan.(ws̄θ̄dz)] .= 0 
    ws̄θ̄dz_avg_top = sum(ws̄θ̄dz[lvls[1]-1, :] .* ΔY_mask) ./ sum(ΔY_mask)
    ws̄θ̄dz_avg_bot = sum(ws̄θ̄dz[lvls[end]+1, :] .* ΔY_mask) ./ sum(ΔY_mask)

    push!(Vdθdy, vs̄θ̄dy_avg)
    push!(Wdθdz_top, ws̄θ̄dz_avg_top)
    push!(Wdθdz_bot, ws̄θ̄dz_avg_bot)

end

l = @layout [a ; b c]
plot(tecco, 1e2.*Vdθdy, ylabel = "cm/s")
plot!(tecco, Wdθdz_bot .- Wdθdz_top, ylabel = "cm/s")
savefig(p, plotsdir("Control_volume_velocites"))




dθ_approx = Vdθdy .- Wdθdz_top .+ Wdθdz_bot
# plot(tecco[1:end-1], inv(2) .* (dθ_approx[1:end-1] .+ dθ_approx[2:end] ) )

θ_approx =  cumsum(cat([0.002], dθ_approx.* 2592000, dims = 1))
plot(tecco, -θ_approx[1:end-1])
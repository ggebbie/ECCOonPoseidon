#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LaTeXStrings, PyCall, GibbsSeaWater 
import PyPlot as plt, NaNMath as nm
@pyimport seaborn as sns
@pyimport pandas as pd
# colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
colors =  sns.color_palette("deep")[1:4]
include(srcdir("config_exp.jl"))
cm = pyimport("cmocean.cm");colorway = cm.balance;
sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0,
              palette = sns.color_palette("colorblind"));
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "PAC56"
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "full")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

#computed on standard pressures, 
#may have a different value if using pressures diagnosed from ECCO
#https://ecco-v4-python-tutorial.readthedocs.io/VectorCalculus_ECCO_barotropicVorticity.html
runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)
nz = length(z)
nt = length(tecco)
P = MeshArray(γ,Float32, nz)
for ijk in eachindex(P)
    P[ijk] .= -pstdz[ijk[2]]
end
Y = abs.((z))
X = zonal_sum(ϕ .* area) ./ zonal_sum(area)
X_mask = X[X .> -40]

colors = [cm.delta, cm.dense, cm.balance]
titles = ["Salinity", "Density", "Temperature"]
p₀ = 2000
ρθSavg = Dict()
for expname in ["iter129_bulkformula", "iter0_bulkformula"]
    filelist = searchdir(diagpath[expname],"state_3d_set1") # first filter for state_3d_set1
    datafilelist_θ  = filter(x -> occursin("data",x),filelist) # second filter for "data"

    zonal_avg = Dict()

    start = 36
    ntt = length(tecco) - start + 1

    θz_mean = MeshArray(γ,Float32, nz); fill!(θz_mean, 0.0)
    σz_mean = MeshArray(γ,Float32, nz); fill!(σz_mean, 0.0)
    Sz_mean = MeshArray(γ,Float32, nz); fill!(Sz_mean, 0.0)

    for tt in 36:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        θname = datafilelist_θ[tt]

        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time θSz = γ.read(diagpath[expname]*θname,MeshArray(γ,Float32,2*nz))
        θz = θSz[:, 1:nz]; Sz = θSz[:, nz+1:end]
        for ijk in eachindex(P)
            σ = OHC_helper.densityJMD95.(θz.f[ijk],Sz.f[ijk], P[ijk], p₀) #EOS from MITGCM 
            σ = σ .- 1000
            θz_mean[ijk] .+= θz.f[ijk] ./ ntt
            Sz_mean[ijk] .+= Sz.f[ijk] ./ ntt
            σz_mean[ijk] .+= σ ./ ntt
        end 

    end

    θ_zonal = OHC_helper.ma_zonal_avg(θz_mean, cell_volumes)
    sigma_zonal = OHC_helper.ma_zonal_avg(σz_mean, cell_volumes)
    S_zonal = OHC_helper.ma_zonal_avg(Sz_mean, cell_volumes)

    zonal_avg["θ"] =  θ_zonal
    zonal_avg["σ2"] =  sigma_zonal
    zonal_avg["S"] =  S_zonal

    ρθSavg[expname] = zonal_avg
end

svename = datadir("native/native_sigma2_zonalavg_" * region *"_1995_2017.jld2")
jldsave(svename, ρθSavg = ρθSavg)

ρθSavg = jldopen(svename)["ρθSavg"]

fig,axes=plt.subplots(1,1, sharey = true, figsize = (15, 20))
time_mean_zonal_avg = ρθSavg["iter129_bulkformula"]["σ2"]
bounds = nm.extrema(time_mean_zonal_avg)
println(bounds)
CS = axes.contourf(X_mask, Y, time_mean_zonal_avg[:, X .> -40], 
cmap=colors[2], vmin = minimum(ECCOtour.sigma2grid()), vmax = maximum(ECCOtour.sigma2grid()), 
levels = ECCOtour.sigma2grid(), extend = "both");
CS = axes.contour(X_mask, Y, time_mean_zonal_avg[:, X .> -40], colors="black", levels = ECCOtour.sigma2grid());
axes.clabel(CS, fontsize=20, inline=true)
axes.set_xlim(-40, 60)
# axes.set_title("Density in ECCO (" * expname* ")")
axes.invert_yaxis()
axes.set_xticks(-40:10:60)
axes.set_xlim(-39, 60)
fig
fig.savefig(plotsdir("NativePlots/DensityAvg" * "_" * expname * ".png"), dpi = 200)
# end

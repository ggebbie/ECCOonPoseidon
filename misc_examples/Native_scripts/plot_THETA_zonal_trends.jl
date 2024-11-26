#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots, PyCall
import NaNMath as nm
import PyPlot as plt
include(srcdir("config_exp.jl"))
cm = pyimport("cmocean.cm");
@pyimport seaborn as sns;

sns.set_theme(context = "poster", style = "ticks", font_scale = 1.0, palette = sns.color_palette("colorblind"));
(ϕ,λ) = latlonC(γ)
area = readarea(γ)

region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
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

nz = length(z[:])
nt = length(tecco)

X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)
Y = -z[:]

for expname in keys(shortnames)
    # Get list of files for salinity on sigma1
    filelist = searchdir(diagpath[expname],"state_3d_set1") 
    datafilelist_θ  = filter(x -> occursin("data",x),filelist)

    println("\n obtaining timeseries for experiment " * expname * " in " * region * "\n")
    θ_zonal_anom_trends = zeros(nz, 270)

    @time for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        fnameθ = datafilelist_θ[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time θ_gcm = γ.read(diagpath[expname]*fnameθ,MeshArray(γ,Float32,50))
        θ_zonal_avg = OHC_helper.ma_zonal_avg(θ_gcm, cell_volumes)
        θ_zonal_anom_trends .+= 100 .* F[2,tt] .* θ_zonal_avg
    end

    boundz = 0.2
    levels = -boundz:0.025:boundz

    fig, ax = plt.subplots(1,1, figsize = (20, 12.5))
    CS = ax.contourf(X, Y, θ_zonal_anom_trends, vmin = -boundz, vmax = boundz, levels = levels, 
    extend = "both", cmap = cm.balance);
    ax.invert_yaxis()
    ax.set_xlabel("Latitude")
    ax.set_ylabel("Depth [m]")
    ax.set_xticks(-56:10:60)
    ax.set_xlim(-49, 60)
    fig.colorbar(CS, orientation = "horizontal", fraction = 0.05, 
    label = L"[°C / century ]")
    fig.savefig(plotsdir("THETA_zonal_trends_" * expname* ".png"), bbox_inches="tight",
    dpi = 1000)
end
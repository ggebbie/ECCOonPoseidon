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
include("./IsopycnalHelpers.jl")
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
tstart = 12*3; tstop = 312; nt = tstop - tstart + 1;
E,F = trend_matrices(tecco[tstart:tstop])

sig1grid = sigma2grid()
nσ = length(sig1grid)
nt = length(tecco)
TSroot = "THETA_on_sigma2" 

[PAC_msk[ff][PAC_msk[ff] .== 0.0] .= NaN for ff = 1:5]

# for expname in keys(shortnames)
expname = "iter129_bulkformula"
# Get list of files for salinity on sigma1
filelist = searchdir(runpath[expname]*"sigma2/",TSroot) # first filter for state_3d_set1
datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
println("\n obtaining timeseries for experiment " * expname * " in " * region * "\n")

# θ_zonal = zeros(nσ, 270, nt)
θ_zonal_anom_trends = zeros(nσ, 270)
@time for (i, tt) in enumerate(tstart:tstop)
    println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
    Tname = datafilelist[tt]
    # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
    @time θσ2 = γ.read(runpath[expname]*"sigma2/"*Tname,MeshArray(γ,Float32,nσ))
    θ_zonal_anom_trends .+= 100 .* F[2,i] .* zonal_avg(θσ2 .* PAC_msk) #C/century
end

# fname = expname * region * "_AVG_P_sigma1.jld2" 
# Pσ = load(σ1datadir(fname), "P")
# sigma2p = mean(Pσ, dims = 2)

# p2z = gsw_z_from_p.(abs.(sigma2p), 30., 0., 0.)
# θ_zonal_anom_trends =  θ_zonal_anom_trends[(!isnan).(p2z[:, 1]), :]
# Y = sig1grid[(!isnan).(p2z[:, 1])]

# p2z = p2z[(!isnan).(p2z[:, 1])]

# iσ = findall(-3500 .< p2z .< -400); 

# θ_zonal_anom_trends = θ_zonal_anom_trends[iσ, :]

Y = sigma2grid()
X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)


boundz = 0.2
levels = -boundz:0.025:boundz

fig, ax = plt.subplots(1,1, figsize = (20, 12.5))
CS = ax.contourf(X, Y, θ_zonal_anom_trends, vmin = -boundz, vmax = boundz, levels = levels, 
extend = "both", cmap = cm.balance);
ax.invert_yaxis()
ax.set_xlabel("Latitude")
ax.set_ylabel("σ₂")
ax.set_xticks(-56:10:60)
ax.set_xlim(-49, 60)
fig.colorbar(CS, orientation = "horizontal", fraction = 0.05, 
label = L"[°C / century ]")
fig
    # fig.savefig(plotsdir("sigma_watermass_trends_" * expname* ".png"), bbox_inches="tight",
    # dpi = 1000)
# endå
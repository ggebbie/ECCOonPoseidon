#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, JLD2,
NCDatasets, Printf, 
DataFrames, LaTeXStrings,
Plots, PyCall, Interpolations
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
E,F = trend_matrices(tecco)

sig1grid = sigma1grid()
nσ = length(sig1grid)
nt = length(tecco)
TSroot = "THETA_on_sigma1" 

for ff = 1:5
    PAC_msk[ff][PAC_msk[ff] .== 0.0] .= NaN
end

Y = sig1grid[40:end]
X = OHC_helper.ma_zonal_sum(ϕ .* area) ./ OHC_helper.ma_zonal_sum(area)
expname = "iter129_bulkformula"
# for expname in keys(shortnames)
    # Get list of files for salinity on sigma1
    filelist = searchdir(runpath[expname]*"sigma1/",TSroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    filelist = searchdir(runpath[expname]*"sigma1/","p_on_sigma1" ) # first filter for state_3d_set1
    datafilelistP  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    println("\n obtaining timeseries for experiment " * expname * " in " * region * "\n")

    # θ_zonal = zeros(nσ, 270, nt)
    θ_zonal_anom_trends = zeros(nσ, 270)
    P_zonal = zeros(nσ, 270)

    @time for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = datafilelist[tt]
        Pname = datafilelistP[tt]

        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))
        θ_zonal_anom_trends .+= 100 .* F[2,tt] .* zonal_avg(θσ1 .* PAC_msk) #C/century
        @time Pσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))
        P_zonal .+=  F[1,tt] .* zonal_avg(Pσ1 .* PAC_msk) #C/century

    end

    
    sigma2p = mean(P_zonal, dims = 2)
    
    p2z = gsw_z_from_p.(abs.(sigma2p), 30., 0., 0.)
    θ_zonal_anom_trends =  θ_zonal_anom_trends[(!isnan).(p2z[:, 1]), :]
    Y = sig1grid[(!isnan).(p2z[:, 1])]
    P_zonal = P_zonal[(!isnan).(p2z[:, 1]), :]

    iz = findall(-3500 .< z .< -400);
    nL = length(X)
    θσ_itp = zeros(length(iz), nL)

    p2z = gsw_z_from_p.(abs.(P_zonal), 30., 0., 0.)

    #take βσ from sigma to depth 
    for il = 1:nL
        itp = interpolate((-p2z[:, il],), θ_zonal_anom_trends[:, il], Gridded(Linear()))
        θσ_itp[:, il] = itp(-vec(z[iz]))
    end

    bounds = 0.2
    levels = -bounds:0.025:bounds

    fig, ax = plt.subplots(1,1, figsize = (20, 12.5))
    CS = ax.contourf(X, -vec(z[iz]), θσ_itp, vmin = -bounds, vmax = bounds, levels = levels, 
    extend = "both", cmap = cm.balance);
    ax.invert_yaxis()
    ax.set_xlabel("Latitude")
    ax.set_ylabel("Depth [m]")

    ax.set_xticks(-56:10:60)
    ax.set_xlim(-49, 60)
    fig.colorbar(CS, orientation = "horizontal", fraction = 0.05, 
    label = L"[°C / century ]")
    fig.savefig(plotsdir("sigma_watermass_trends_onZ_" * expname* ".png"), bbox_inches="tight", dpi = 1000)
# end
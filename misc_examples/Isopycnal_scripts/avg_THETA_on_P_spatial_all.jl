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
import PyPlot as plt
using PyCall
cm = pyimport("cmocean.cm");colorway = cm.balance;

include(srcdir("config_exp.jl"))
include("./IsopycnalHelpers.jl")

#plotting stuff 
proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
cf = Vector{Any}(undef ,1)

#reshape λ for plotting 
(ϕ,λ) = latlonC(γ)
for ff ∈ [3,4] 
    λ[ff][λ[ff] .<= 0 ] = λ[ff][λ[ff] .<= 0 ] .+ 360
end

#directory stuff 
σ1datadir(x="") = OHC_helper.σ1datadir(x)
!ispath(σ1datadir()) && mkpath(σ1datadir())

area = readarea(γ)

region = "PAC56"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
cell_depths = OHC_helper.get_cell_depths(PAC_msk, ΔzF, Γ.hFacC); 
cell_volumes = OHC_helper.get_cell_volumes(area, cell_depths);

runpath,diagpath = listexperiments(exprootdir());
ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)
experiments = keys(shortnames)

tecco= 1992+1/24:1/12:2018 # ecco years
E,F = trend_matrices(tecco)

expname = "iter129_bulkformula"

sig1grid = sigma1grid()
nσ = length(sig1grid)
#get root for theta on sigma1
TSroot = "THETA_on_sigma1" 

β_dict = Dict()
global expname

for (i, expr) in enumerate(experiments)
    expname = expr; println(expname)
    filelist = searchdir(runpath[expname]*"sigma1/",TSroot) # first filter for state_3d_set1
    datafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    nt = length(datafilelist)

    fname = expname * region * "_AVG_P_sigma1.jld2" 
    !ispath(σ1datadir(fname)) && (include("get_avgP_on_sigma1.jl"))
    Pσ = load(σ1datadir(fname), "P")

    sigma2p = mean(Pσ, dims = 2)
    p2z = gsw_z_from_p.(abs.(sigma2p), 30., 0., 0.); p2z = vec(p2z)
    σlvl = findmin(p2z .- 1200.)[2]
    β = MeshArray(γ,Float32); fill!(β, 0.0)

    @time for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = datafilelist[tt]
        # get S on sigma1. Way to read a slice? (Didn't figure it out yet)
        @time θσ1 = γ.read(runpath[expname]*"sigma1/"*Tname,MeshArray(γ,Float32,nσ))
        θσ1 = θσ1[:, σlvl]
        for ff = 1:5
            β.f[ff] .+= F[2,tt] .* θσ1.f[ff] 
        end
    end
    β_dict[expr] = β
end

#plotting climatology :-)
fig, axes = plt.subplots(2, 2, figsize=(15,10), subplot_kw=Dict("projection"=> proj0))
# bounds = maximum([nm.maximum([nm.maximum(β[ff].*PAC_msk[ff]) for ff = 1:5]) for β in values(β_dict)])
labels = ["Forcing + Init. Cond. ", "Forcing", "Init. Cond.", "Neither"]

for (i, expr) in enumerate(experiments)
    ax = axes[i]
    β = β_dict[expr]
    println(labels[i])
    β_mask = deepcopy(β); β_mask[findall(β_mask.==0)] = NaN
    ax.coastlines(resolution="110m")
    ax.set_extent((-250, -80, -40, 70),crs=projPC)
    gl = ax.gridlines(crs=projPC, draw_labels=true,
                        linewidth=2, color="gray", alpha=0, linestyle="--")
    gl.xlabels_top = false;gl.ylabels_right = false

    for ff in 1:5
        cf[1] = ax.pcolormesh(λ[ff], ϕ[ff],  β_mask[ff],
        vmin = -bounds, vmax = bounds, shading="nearest", transform=projPC, 
        rasterized = true, cmap = colorway)               
    end
    # ax.set_title(labels[i])
    ax.tick_params(bottom=true, left=true)
    ax.set_title(labels[i])

end
fig.tight_layout()
# fig.suptitle("Temperature Trends between 2-3 km depth on sigma1)")
fig
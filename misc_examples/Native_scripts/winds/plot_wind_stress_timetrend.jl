#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../src/intro.jl")
include("../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall
import PyPlot as plt
import NaNMath as nm
@pyimport matplotlib.animation as anim

cm = pyimport("cmocean.cm");
@pyimport seaborn as sns;

sns.set_theme(context = "paper", style = "ticks",
              palette = sns.color_palette("colorblind"));
              
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)
area = readarea(γ)

ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")

tecco = 1992+1/24:1/12:2018
nz = 50; nt = 312
E,F = trend_matrices(tecco)

runpath,diagpath = listexperiments(exprootdir());
τx_dict = Dict(); τy_dict = Dict()
vars = ["iter129_bulkformula", "seasonalclimatology"]
for expname in vars
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    τx_dict[expname] = MeshArray(γ,Float32); fill!(τx_dict[expname], 0.0)
    τy_dict[expname] = MeshArray(γ,Float32); fill!(τy_dict[expname], 0.0)
    @time for (i, tt) = enumerate(1:nt)
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = τdatafilelist[tt]
        
        τx, τy = OHC_helper.extract_ocnTAU(diagpath, expname , τdatafilelist, tt, γ)
        τxC, τyC = velocity2centerfast(τx, τy, Γ)
        τE, τN = rotate_uv(τxC, τyC, Γ)

        for ff = 1:5
            τx_dict[expname].f[ff] .+= F[2,tt] .* τE.f[ff] .* 100 # per century
            τy_dict[expname].f[ff] .+= F[2,tt] .* τN.f[ff] .* 100 # per century
        end
    end
end

diff_dict(d, x, y) = d[x] .- d[y]
diff_label = shortnames[vars[1]] * " minus " * shortnames[vars[2]]
τx_dict[diff_label] = diff_dict(τx_dict, vars[1], vars[2])
τy_dict[diff_label] = diff_dict(τy_dict, vars[1], vars[2])
push!(vars, diff_label) 

proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
λ_wrap = OHC_helper.wrap_λ(λ)

#setting up plotting stuff       
bounds = [0.15, 0.15, 0.15]
CF = Any[1, 1, 1]

#plot τx
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)

for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, τx_dict[var], 
                        bounds[i], cm.curl, projPC, true)
    ax[i].set_title(var)
end

[a.coastlines() for a in ax]
[a.set_extent((-326, -60, -60, 70),crs=projPC) for a in ax]
fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[N/m²] per century")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[N/m²] per century")
fig.suptitle("Zonal Windstress Trend[1993 - 2017]", y = 0.25)
fig
fig.savefig(plotsdir("native/trendτx_" * shortnames[vars[1]] * "_" * shortnames[vars[2]] * ".png"), bbox_inches = "tight", dpi = 400);

#plot τy
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
                       
for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, τy_dict[var], 
                        bounds[i], cm.curl, projPC, true)
    ax[i].set_title(var)
end

[a.coastlines() for a in ax]
[a.set_extent((-326, -60, -60, 70),crs=projPC) for a in ax]
fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[N/m²] per century")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[N/m²] per century")
fig.suptitle("Meridional Windstress Trend[1993 - 2017]", y = 0.25)
fig.savefig(plotsdir("native/trendτy_" * shortnames[vars[1]] * "_" * shortnames[vars[2]] * ".png"), bbox_inches = "tight", dpi = 400);
#analysis should complete within 4 minutes 
#using 1 threads 
# julia --threads=4 --project=@. ./extract_theta_native.jl

include("../../src/intro.jl")
include("../../src/OHC_helper.jl")

using Revise,ECCOonPoseidon, ECCOtour,
MeshArrays, MITgcmTools, JLD2, DrWatson, Statistics, 
Printf, LinearAlgebra, LaTeXStrings, PyCall
import PyPlot as plt
import NaNMath as nm
@pyimport matplotlib.animation as anim

cm = pyimport("cmocean.cm");
@pyimport seaborn as sns;

sns.set_theme(context = "notebook", style = "ticks",
              palette = sns.color_palette("colorblind"));
              
include(srcdir("config_exp.jl"))

(ϕ,λ) = latlonC(γ)

ignore_list= ["noIA", "129ff"]
shortnames = OHC_helper.reduce_dict(expnames(), ignore_list)

region = "PAC"; 
PAC_msk = OHC_helper.PAC_mask(Γ, basins, basin_list, ϕ, λ; 
region, extent = "not")
ocean_mask = OHC_helper.wet_pts(Γ)

tecco= 1992+1/24:1/12:2018 # ecco years
nt = length(tecco)

runpath,diagpath = listexperiments(exprootdir());
diagpath["clim_tau_iter0"] = "/vast/ECCOv4r4/exps/clim_tau_iter0/run/diags/"

τx_dict = Dict(); τy_dict = Dict(); curlτ_dict = Dict(); ekup_dict = Dict()
vars = ["iter0_bulkformula", "iter129_bulkformula", "only_init", "only_kappa", "only_sfc"]

for expname in vars
    filelist = searchdir(diagpath[expname],"state_2d_set1") # first filter for state_3d_set1
    τdatafilelist  = filter(x -> occursin("data",x),filelist) # second filter for "data"
    τx_dict[expname] = MeshArray(γ,Float32); fill!(τx_dict[expname], 0.0)
    τy_dict[expname] = MeshArray(γ,Float32); fill!(τy_dict[expname], 0.0)

    for tt = 1:nt
        println("year ",Int(floor(tecco[tt]))," month ",((tt-1)%12)+1)
        Tname = τdatafilelist[tt]
        
        τx, τy = OHC_helper.extract_ocnTAU(diagpath, expname , τdatafilelist, tt, γ)

        for ff = 1:5
            τx_dict[expname].f[ff] .+= τx.f[ff] ./ nt
            τy_dict[expname].f[ff] .+= τy.f[ff] ./ nt
        end
    end
    curlτ_dict[expname] = MeshArrays.curl(τx_dict[expname], τy_dict[expname], Γ)
    τxC, τyC = velocity2center(τx_dict[expname], τy_dict[expname], Γ)
    τE, τN = rotate_uv(τxC, τyC, Γ)
    τx_dict[expname] = τE
    τy_dict[expname] = τN
end

diff_dict(d, x, y) = d[x] .- d[y]
shortnames["climatological_tau"] = "mean_wind_i129"
shortnames["clim_tau_iter0"] = "mean_wind_i0"

diff_label = shortnames[vars[1]] * " minus " * shortnames[vars[2]]
curlτ_dict[diff_label] = diff_dict(curlτ_dict, vars[1], vars[2])
τx_dict[diff_label] = diff_dict(τx_dict, vars[1], vars[2])
τy_dict[diff_label] = diff_dict(τy_dict, vars[1], vars[2])

push!(vars, diff_label) 

# proj0 = ECCOonPoseidon.cartopy.crs.PlateCarree(central_longitude=-150)
proj0 = ECCOonPoseidon.cartopy.crs.Robinson(central_longitude=-150)
projPC = ECCOonPoseidon.cartopy.crs.PlateCarree()
λ_wrap = OHC_helper.wrap_λ(λ)

#setting up plotting stuff       
bounds = [0.15, 0.15, 0.04]
CF = Any[1, 1, 1]

#plot τx
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)

for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, τx_dict[var], 
                        bounds[i], cm.curl, projPC, true)
    ax[i].set_title(var)
end
# ax.set_extent((-180, 180, -70, 56),crs=projPC)

[a.coastlines() for a in ax]
[a.set_extent((-180, 180, -70, 56),crs=projPC) for a in ax]
fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[N/m²]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[N/m²]")
fig.suptitle("Mean Zonal Wind Stress [1993 - 2017]", y = 0.25)
fig
# fig.savefig(plotsdir("native/mean_τx_" * shortnames[vars[1]] * "_" * shortnames[vars[2]] * ".png"), bbox_inches = "tight", dpi = 400);

#plot τy
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
                       
for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, τy_dict[var], 
                        bounds[i], cm.curl, projPC, true)
    ax[i].set_title(var)
end

[a.coastlines() for a in ax]
[a.set_extent((-180, 180, -70, 56),crs=projPC) for a in ax]
fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[N/m²]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[N/m²]")
fig.suptitle("Mean Meridional Wind Stress [1993 - 2017]", y = 0.25)
fig
# fig.savefig(plotsdir("native/mean_τy_" * shortnames[vars[1]] * "_" * shortnames[vars[2]] * ".png"), bbox_inches = "tight", dpi = 400);


#plot curl
fig, ax = plt.subplots(1, 3, figsize=(15,20), 
                       subplot_kw=Dict("projection"=> proj0), sharey = true)
[a.coastlines() for a in ax]
[a.set_extent((120, 285, -70, 56),crs=projPC) for a in ax]

for (i, var) in enumerate(vars)
    CF[i] = OHC_helper.pcolormesh_ma(ax[i], λ_wrap, ϕ, curlτ_dict[var], 
                        1e-12, cm.curl, projPC, true)
    ax[i].set_title(var)
end

fig.colorbar(CF[1], ax = ax[1:2], fraction = 0.01, pad = 0.03, orientation = "horizontal", label = "[N/m³]")
fig.colorbar(CF[3], ax = ax[3], fraction = 0.01,  pad = 0.03, orientation = "horizontal", label = "[N/m³]")
fig.suptitle("Mean Wind Stress Curl [1992 - 2017]", y = 0.31)
# fig.savefig(plotsdir("native/mean_curlτ_" * shortnames[vars[1]] * "_" * shortnames[vars[2]] * ".png"), bbox_inches = "tight", dpi = 400);
fig